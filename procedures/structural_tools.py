from parsers.icsv import ICSV
import logging
import logging.handlers
import procedures.kidera
import procedures.pdbencode
import procedures.popscomp
import procedures.contacts as contacts
import json
from Bio.SeqUtils import seq1
import string
from Bio.PDB import PDBParser
from Bio import BiopythonWarning
import warnings
import gzip
import shutil
import os
import procedures.common as common
import pickle as pickle
from concurrent.futures import ProcessPoolExecutor
import multiprocessing as multi
import concurrent.futures as cf
import threading

# Multiprocess implementation with single logging sink
# Forced spawn method of process creation due to issues with common logging
#
# References:
#   https://docs.python.org/dev/howto/logging-cookbook.html#logging-to-a-single-file-from-multiple-processes
#   https://docs.python.org/3/library/multiprocessing.html#contexts-and-start-methods
#   https://superfastpython.com/processpoolexecutor-in-python/
#   https://rpy2.github.io/doc/latest/html/callbacks.html
#
# Note: moving to spawning multiprocesses resulted in more memory usage in Linux
#   so 4GB of RAM memory is recommended


data_results_filename_suffix = ".data_results.pickle"

# notice some overlaps - definitions from Ittisoponpisan S et al.
res_hydrophobic = ["ALA", "VAL", "ILE", "GLY", "LEU", "MET"]
res_hydrophilic = ["SER", "THR", "CYS", "PRO", "GLN", "ASN"]
res_negative = ["ASP", "GLU"]
res_positive = ["ARG", "HIS", "LYS"]
res_aromatic = ["TYR", "TRP", "PHE"]


# @dataclass
# class TaskResult:
#     popscomp_res: bool = False
#     popscomp_err: str = ""
#     general_res: bool = False
#     general_err: str = ""


def _make_task_result() -> dict:
    return {
        "popscomp_res": False,
        "popscomp_err": "",
        "processentry_res": False,
        "processentry_err": "",
        "contacts_list_res": False,
        "contacts_list_err": "",
        "general_res": False,
        "general_err": "",
    }


def logging_thread(queue: multi.Queue):
    while True:
        record = queue.get()
        if record is None:
            break
        logger = logging.getLogger()
        logger.handle(record)


def subprocess_logging_configurer(queue: multi.Queue):
    h = logging.handlers.QueueHandler(queue)
    root = logging.getLogger()
    root.addHandler(h)
    root.setLevel(logging.DEBUG)


class StructuralToolsRunner:
    def __init__(
        self,
        node_pdb_filename: str,
        pdb_folder: str,
        structure_db_type: str,
        scrap_folder: str,
        issues_filename: str,
        pdb_residues_progress_filename: str,
        skip_defects: bool,
        homodimers: bool,
        use_cache: bool,
        num_threads=0,
    ):
        self.node_pdb_filename = node_pdb_filename
        self.pdb_folder = pdb_folder
        self.structure_db_type = structure_db_type
        self.scrap_folder = scrap_folder
        self.issues_filename = issues_filename
        self.pdb_residues_progress_filename = pdb_residues_progress_filename
        self.homodimers = homodimers
        self.items_read = 0
        self.num_items = 0
        self.skip_defects = skip_defects
        self.use_cache = use_cache
        self.num_threads = 0
        self.new_progress_set = set()
        if num_threads < 0:
            self.num_threads = 1
        elif num_threads == 0:
            self.num_threads = multi.cpu_count()
        elif num_threads > multi.cpu_count():
            self.num_threads = multi.cpu_count()
        else:
            self.num_threads = num_threads

    def _load_issues_from_file(self):
        self.issues = dict()
        self.issues["items_read"] = 0
        self.issues["max_items"] = -1
        self.issues["report"] = dict()
        # read issues file
        try:
            issues_file = open(self.issues_filename, mode="r")
            self.issues = json.load(issues_file)
            issues_file.close()
        except json.JSONDecodeError as j:
            logging.fatal("Error whilst parsing existing issue file: %s", j)
            exit(code=1)
        except Exception:
            logging.info("Issues file not found, creating a new one")

        self.skip_items = self.issues["items_read"]
        self.max_items = self.issues["max_items"]

    def _store_run_info(self):
        self.issues["items_read"] = self.items_read
        # num_def_physio = 0
        # num_def_nonphysio = 0
        # for code in self.issues["defects"]:
        #     if self.issues["defects"][code]["physio"]:
        #         num_def_physio = num_def_physio + 1
        #     else:
        #         num_def_nonphysio = num_def_nonphysio + 1
        # self.issues["num_def_physio"] = num_def_physio
        # self.issues["num_def_nonphysio"] = num_def_nonphysio
        # self.issues["num_def"] = num_def_nonphysio + num_def_physio
        issues_file = open(self.issues_filename, mode="w")
        json.dump(self.issues, issues_file, indent=True)
        issues_file.close()

        # self.pdb_residues_progress_file.flush()
        (
            pdb_residues_progress_file,
            pdb_residues_progress,
        ) = common.create_or_append_file(
            self.pdb_residues_progress_filename, [self.progress_dict_key]
        )
        for item in self.new_progress_set:
            pdb_residues_progress.writerow([item])
        pdb_residues_progress_file.close()
        self.new_progress_set.clear()

        # # export DBs to JSON
        # physio_file = open(self.physio_filename, mode="w")
        # json.dump(self.db_physio, physio_file, indent=True)
        # physio_file.close()

        # nonphysio_file = open(self.nonphysio_filename, mode="w")
        # json.dump(self.db_nonphysio, nonphysio_file, indent=True)
        # nonphysio_file.close()

    def _select_dispatch_items(self):
        logging.info(
            "Starting running select and dispatch items with {} threads".format(
                self.num_threads
            )
        )

        queue = multi.Manager().Queue(-1)
        lp = threading.Thread(target=logging_thread, args=(queue,))
        # start logger
        lp.start()

        # force spawn of processes rather than fork as in Linux
        # forking causes issues with common logging
        spawn_context = multi.get_context("spawn")

        with ProcessPoolExecutor(
            max_workers=self.num_threads,
            mp_context=spawn_context,
            initializer=subprocess_logging_configurer,
            initargs=(queue,),
        ) as executor:
            try:
                # submit tasks and collect futures
                futures = [
                    executor.submit(self._run_single_item, pdb_code)
                    for pdb_code in self.node_pdb
                ]

                # process task results as they are available
                for future in cf.as_completed(futures):
                    self.items_read += 1
                    # retrieve the result
                    result = future.result()
                    result_pdb_code = result["pdb_code"]
                    logging.info(
                        f"Processed {result_pdb_code} ({self.items_read}/{self.num_items})"
                    )
                    if result["task_result"] is not None:
                        if not result["task_result"]["general_res"]:
                            self.issues["report"][result_pdb_code] = result[
                                "task_result"
                            ]

                        # processed this PDB so save progress
                        self.new_progress_set.add(result_pdb_code)

                        if self.items_read % 10 == 0:
                            self._store_run_info()
            except Exception as e:
                logging.error(f'Got exception "{e}", waiting')
                executor.shutdown(wait=True, cancel_futures=True)
                logging.error("Finished waiting")
                queue.put(None)
                lp.join()

        queue.put(None)
        lp.join()

        return

    # this method runs in different processes from ProcessPoolExecutor
    # make sure that shared data is not written here
    def _run_single_item(self, pdb_code: str):
        # basic checks for skipping items
        if pdb_code in self.previous_progress:
            logging.info("Skipping previously processed entry: {}".format(pdb_code))
            return {"task_result": None, "pdb_code": pdb_code}

        if self.skip_defects and pdb_code in self.issues["report"]:
            if self.issues["report"][pdb_code]["general_res"] is False:
                logging.info("Skipping entry with known issue: {}".format(pdb_code))
                return {"task_result": None, "pdb_code": pdb_code}

        logging.info(f"Processing {pdb_code}")
        task = StructuralToolsTask(
            pdb_code=pdb_code,
            pdb_folder=self.pdb_folder,
            structure_db_type=self.structure_db_type,
            scrap_folder=self.scrap_folder,
            homodimers=self.homodimers,
            use_cache=self.use_cache,
        )
        result = {"task_result": task.run(), "pdb_code": pdb_code}
        return result

    def run(self):
        node_file_type = ""
        node_dict_key = ""
        progress_file_type = ""
        self.progress_dict_key = ""
        match self.structure_db_type:
            case "pdb" | "biounit":
                node_file_type = "node_pdb"
                node_dict_key = "id:ID(pdb-id)"
                progress_file_type = "pdb_residues_progress"
                self.progress_dict_key = "pdb-id"

            case "alphafold":
                node_file_type = "node_model"
                node_dict_key = "id:ID(model-id)"
                progress_file_type = "model_residues_progress"
                self.progress_dict_key = "model-id"

        self.node_pdb = _load_dictionary_from_csv_file(
            file_name=self.node_pdb_filename,
            file_type=node_file_type,
            dict_key=node_dict_key,
        )
        progress_dict = _load_dictionary_from_csv_file(
            file_name=self.pdb_residues_progress_filename,
            file_type=progress_file_type,
            dict_key=self.progress_dict_key,
            file_not_found_fails=False,
        )
        self.previous_progress = set(progress_dict.keys())

        self.num_items = len(self.node_pdb)

        self._load_issues_from_file()

        self._select_dispatch_items()

        self._store_run_info()


class StructuralToolsTask:
    def __init__(
        self,
        pdb_code: str,
        pdb_folder: str,
        structure_db_type: str,
        scrap_folder: str,
        homodimers: bool,
        use_cache: bool,
    ):
        self.pdb_code = pdb_code
        self.subfolder_divider_code = self.pdb_code[1:3]
        self.pdb_folder = pdb_folder
        self.structure_db_type = structure_db_type
        self.scrap_folder = scrap_folder
        self.pdb_file = ""
        self.scrap_pdb_folder = "{0}/{1}/{2}".format(
            self.scrap_folder, self.subfolder_divider_code, self.pdb_code
        )
        self.homodimers = homodimers
        self.data_results_filename = "{0}/{1}{2}".format(
            self.scrap_pdb_folder, self.pdb_code, data_results_filename_suffix
        )

        # self.max_pdb_file_size = 256000  # 250 KB
        self.max_pdb_file_size = 1 * 1048576  # 10 MB
        self.popscomp = procedures.popscomp.POPScomp(
            scrap_folder=self.scrap_folder, visualise=False
        )
        self.pdbe = procedures.pdbencode.PDBEncode("")
        self.pops_residues = dict()
        self.region_data = dict()
        self.contacts_list = list()
        # TODO: remove this
        self.is_physio = False
        self.db_nonphysio = {
            "dimers": dict(),
        }
        self.db_physio = {
            "dimers": dict(),
        }
        self.issues = dict()
        self.use_cache = use_cache
        self.result = _make_task_result()

    def run(self) -> dict:
        # seqfull = item["seqfull"]
        # len_seq = len(seqfull)

        # MOVE OUT
        # self._set_output_files()

        if self.use_cache:
            self._run_cached()
        else:
            self._run_no_cached()

        self.result["general_res"] = (
            self.result["popscomp_res"]
            and self.result["processentry_res"]
            and self.result["contacts_list_res"]
        )

        self._write_data_results_file()

        # MOVE OUT
        # # output to csv files
        # self._output_chain_residue_to_csv(region_data, self.pdb_code)

        self._clean_up(self.pdb_file, self.scrap_pdb_folder)

        # MOVE OUT
        # self._close()
        return self.result

    def _run_cached(self):
        self._load_data_results_file()
        if "general_res" not in self.result or not self.result["general_res"]:
            if not self._prepare_pdb_file():
                return

        if "popscomp_res" not in self.result or not self.result["popscomp_res"]:
            if not self._run_popscomp():
                return

        if "processentry_res" not in self.result or not self.result["processentry_res"]:
            if not self._run_process_entry():
                return
            if self.structure_db_type == "alphafold":
                if not self._run_get_alphafold_pLDDT():
                    return

        if (
            "contacts_list_res" not in self.result
            or not self.result["contacts_list_res"]
        ):
            if not self._run_contact_list():
                return

    def _run_no_cached(self):
        if not self._prepare_pdb_file():
            return

        if not self._run_popscomp():
            return

        if not self._run_process_entry():
            return

        if self.structure_db_type == "alphafold":
            if not self._run_get_alphafold_pLDDT():
                return

        if not self._run_contact_list():
            return

    # def _set_output_files(self):
    #     # (
    #     #     self.pdb_residues_progress_file,
    #     #     self.pdb_residues_progress,
    #     # ) = common.create_or_append_file(
    #     #     self.pdb_residues_progress_filename, ["pdb-id"]
    #     # )
    #     self.node_chain_file, self.node_chain_csv = common.create_or_append_file(
    #         self.node_chain_filename,
    #         [
    #             "id:ID(chain-id)",
    #             "name:string",
    #             "core_numres:long",
    #             "core_k1:double",
    #             "core_k2:double",
    #             "core_k3:double",
    #             "core_k4:double",
    #             "sface_numres:long",
    #             "sface_k1:double",
    #             "sface_k2:double",
    #             "sface_k3:double",
    #             "sface_k4:double",
    #             "iface_numres:long",
    #             "iface_k1:double",
    #             "iface_k2:double",
    #             "iface_k3:double",
    #             "iface_k4:double",
    #         ],
    #         force_create=True,
    #     )
    #     self.node_residue_file, self.node_residue_csv = common.create_or_append_file(
    #         self.node_residue_filename,
    #         [
    #             "id:ID(residue-id)",
    #             "name:string",
    #             "num:long",
    #             "q_sasa:double",
    #             "location",
    #         ],
    #         force_create=True,
    #     )
    #     self.rel_pdb_chain_file, self.rel_pdb_chain_csv = common.create_or_append_file(
    #         self.rel_pdb_chain_filename,
    #         [
    #             ":START_ID(pdb-id)",
    #             ":END_ID(chain-id)",
    #             ":TYPE",
    #         ],
    #         force_create=True,
    #     )
    #     (
    #         self.rel_pdb_residue_file,
    #         self.rel_pdb_residue_csv,
    #     ) = common.create_or_append_file(
    #         self.rel_pdb_residue_filename,
    #         [
    #             ":START_ID(pdb-id)",
    #             ":END_ID(residue-id)",
    #             ":TYPE",
    #         ],
    #         force_create=True,
    #     )
    #     (
    #         self.rel_chain_residue_file,
    #         self.rel_chain_residue_csv,
    #     ) = common.create_or_append_file(
    #         self.rel_chain_residue_filename,
    #         [
    #             ":START_ID(chain-id)",
    #             ":END_ID(residue-id)",
    #             ":TYPE",
    #         ],
    #         force_create=True,
    #     )

    # def _close(self):
    #     # self.pdb_residues_progress_file.close()
    #     self.node_chain_file.close()
    #     self.node_residue_file.close()
    #     self.rel_pdb_chain_file.close()
    #     self.rel_chain_residue_file.close()
    #     self.rel_pdb_residue_file.close()
    #     # super().close()

    def _prepare_pdb_file(self) -> bool:
        self.pdb_file = ""
        # points to the PDB file copied by popscomp module

        # create entry's scrap directory
        os.makedirs(name=self.scrap_pdb_folder, exist_ok=True)

        if self.homodimers:
            # use the PDB file as it is
            self.pdb_file = self.pdb_folder + "/" + self.pdb_code + ".pdb"
            return True
        else:
            # locate the compressed file, decompress and copy
            # to scrap folder temporarily
            pdb_file_gz = ""
            match self.structure_db_type:
                case "pdb":
                    subfolder = self.pdb_code[1:3]
                    pdb_file_gz = "{}/{}/pdb{}.ent.gz".format(
                        self.pdb_folder, subfolder, self.pdb_code
                    )
                case "biounit":
                    subfolder = self.pdb_code[1:3]
                    pdb_file_gz = "{}/{}/{}.pdb1.gz".format(
                        self.pdb_folder, subfolder, self.pdb_code
                    )
                case "alphafold":
                    pdb_file_gz = "{}/AF-{}-F1-model_v4.pdb.gz".format(
                        self.pdb_folder, self.pdb_code
                    )
                case _:
                    raise Exception(
                        "Unexpected structure db type. Valid values are: pdb, biounit, alphafold"
                    )
            self.pdb_file = self.scrap_folder + "/" + self.pdb_code + ".pdb"
            try:
                # test if file too big and skip
                file_size = os.path.getsize(pdb_file_gz)
                if file_size > self.max_pdb_file_size:
                    raise Exception("PDB file is too big")
                with gzip.open(pdb_file_gz, "rb") as f_in:
                    with open(self.pdb_file, "wb") as f_out:
                        shutil.copyfileobj(f_in, f_out)
                return True
            except Exception as e:
                self.result["general_err"] = "Unable to open PDB file: " + str(e)
                self._clean_up(self.pdb_file, self.scrap_pdb_folder)
                return False

    def _load_data_results_file(self):
        self.pops_residues = dict()
        self.region_data = dict()
        self.result = _make_task_result()
        if os.path.exists(self.data_results_filename):
            with open(self.data_results_filename, "rb") as file:
                data_results = pickle.load(file)
                if "results" in data_results:
                    self.result = data_results["results"]
                if "pops_residues" in data_results:
                    self.pops_residues = data_results["pops_residues"]
                if "region_data" in data_results:
                    self.region_data = data_results["region_data"]
                if "contacts_list" in data_results:
                    self.contacts_list = data_results["contacts_list"]

    def _write_data_results_file(self):
        with open(self.data_results_filename, "wb") as file:
            data_results = {
                "pops_residues": self.pops_residues,
                "region_data": self.region_data,
                "results": self.result,
                "contacts_list": self.contacts_list,
            }
            pickle.dump(obj=data_results, file=file)

    def _run_popscomp(self) -> bool:
        msg = ""
        self.pops_residues, msg = self.popscomp.run(self.pdb_file, self.pdb_code)
        if self.pops_residues is None:  # in fact an exception
            self.result["popscomp_err"] = "Skipped: " + msg
            return False
        if self.homodimers:
            if len(self.pops_residues) != 2:  # other two chains
                self.result["popscomp_err"] = (
                    """Skipped: item does not have exactly
                two chains, structure has {} chains""".format(
                        len(self.pops_residues)
                    )
                )
                return False

        if len(msg) > 0:
            self.result["popscomp_err"] = "Non-fatal: " + msg

        self.result["popscomp_res"] = True
        return True

    def _run_process_entry(self) -> bool:
        # TODO: divide into more steps

        # prepare pdbencode
        self.pdbe.load_structure(file_path=self.pdb_file)

        iface_entry, skip_item, err_msg = self._process_entry(
            "interface", store_residue_data=False
        )
        if skip_item:
            self.result["processentry_err"] = "Skipped item (iface). " + err_msg
            return False
        sface_entry, skip_item, err_msg = self._process_entry(
            "surface", store_residue_data=False
        )
        if skip_item:
            self.result["processentry_err"] = "Skipped item (sface). " + err_msg
            return False
        core_entry, skip_item, err_msg = self._process_entry(
            "core", store_residue_data=False
        )
        if skip_item:
            self.result["processentry_err"] = "Skipped item (core). " + err_msg
            return False

        if self.is_physio:
            # self.db_physio["total_iface_residues"] += entry["num_iface_residues"]
            # self.db_physio["total_hydrophobic"] += entry["num_iface_hydrophobic"]
            # self.db_physio["total_hydrophilic"] += entry["num_iface_hydrophilic"]
            # self.db_physio["total_negative"] += entry["num_iface_negative"]
            # self.db_physio["total_positive"] += entry["num_iface_positive"]
            # self.db_physio["total_aromatic"] += entry["num_iface_aromatic"]
            # self.db_physio["total_unknown"] += entry["num_iface_unknown"]
            # self.db_physio["all_iface_pdbencode"] += entry["iface_pdbencode"]
            self.db_physio["dimers"][self.pdb_code] = {
                "interface": iface_entry,
                "surface": sface_entry,
                "core": core_entry,
            }
        else:
            # self.db_nonphysio["total_iface_residues"] += entry["num_iface_residues"]
            # self.db_nonphysio["total_hydrophobic"] += entry["num_iface_hydrophobic"]
            # self.db_nonphysio["total_hydrophilic"] += entry["num_iface_hydrophilic"]
            # self.db_nonphysio["total_negative"] += entry["num_iface_negative"]
            # self.db_nonphysio["total_positive"] += entry["num_iface_positive"]
            # self.db_nonphysio["total_aromatic"] += entry["num_iface_aromatic"]
            # self.db_nonphysio["total_unknown"] += entry["num_iface_unknown"]
            # self.db_nonphysio["all_iface_pdbencode"] += entry["iface_pdbencode"]
            self.db_nonphysio["dimers"][self.pdb_code] = {
                "interface": iface_entry,
                "surface": sface_entry,
                "core": core_entry,
            }
        # logging.info("Processed " + pdb_code)
        # self.persistence.create_vertex(label="homodimer", props=item)

        self.region_data = {
            "core": core_entry,
            "surface": sface_entry,
            "interface": iface_entry,
        }

        self.result["processentry_res"] = True
        return True

    def _run_contact_list(self) -> bool:
        try:
            self.contacts_list = contacts.generate_contact_residue_list(
                self.pdb_code, self.pdb_file, contact_map_cutoff=7.0
            )
        except Exception as e:
            self.result["contacts_list_err"] = str(e)
            logging.error(f"Exception whilst calculating contact map: {e}")
            return False

        self.result["contacts_list_res"] = True
        return True

    def _process_entry(
        self,
        group_name: str,
        store_residue_data: bool,
    ):  # returns dict, bool, str
        # print(seqfull)
        # process interface residues
        skip_item = False
        err_msg = ""
        num_group_residues = 0
        group_res_type_count = {
            "num_res_hydrophobic": 0,
            "num_res_hydrophilic": 0,
            "num_res_negative": 0,
            "num_res_positive": 0,
            "num_res_aromatic": 0,
            "num_res_unknown": 0,
        }
        group_pdbencode = ""
        chains = dict()
        for chain in self.pops_residues:
            # encode chain using PDBEncode
            # encoded_chain = self.pdbe.encode_chain(chain=chain)
            # print("Encoded:  " + encoded_chain)
            # print("Sequence: " + seqfull)
            kidera_sum = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
            num_chain_residues = 0
            chain_res_type_count = {
                "num_res_hydrophobic": 0,
                "num_res_hydrophilic": 0,
                "num_res_negative": 0,
                "num_res_positive": 0,
                "num_res_aromatic": 0,
                "num_res_unknown": 0,
            }
            chain_pdbencode = ""
            chains[chain] = dict()
            if store_residue_data:
                chains[chain]["residues"] = dict()

            for res_num in self.pops_residues[chain]["residues"]:
                # check if residue is from a desired group name
                if (
                    self.pops_residues[chain]["residues"][res_num]["location"]
                    != group_name
                ):
                    continue

                num_chain_residues = num_chain_residues + 1
                # compare residues from the interface with the sequence (for sanity)

                res_data = self.pops_residues[chain]["residues"][res_num]
                # print("chain {}, resnum {}, interface resname {}, sequence resname {}".format(chain, res_num, res_data["res_name"], seqfull[res_num-1]))
                # if seq1(res_data["res_name"]) != seqfull[res_num-1]:
                #     raise Exception("Mismatch in residues")

                # add kidera factors to all residues in the interface
                res_code = seq1(res_data["res_name"])
                kid_dframe = True
                if len(res_code) == 0:
                    # no residue code mapped into 3 letter (maybe unknown residue?)
                    err_msg = "Unable to calculate Kidera: unknown residue code. Chain: {}, Residue Number: {}, Residue Code: {}\nSkipping structure due to detected defect".format(
                        chain, res_num, res_data["res_name"]
                    )
                    logging.warning(
                        "Skipping item %s due to detected defect", self.pdb_code
                    )
                    skip_item = True
                    # not storing under defects any more
                    # self.issues["defects"][self.pdb_code] = {
                    #     "reason": err_msg,
                    #     "physio": self.is_physio,
                    # }
                    break

                try:
                    kid_dframe = procedures.kidera.score_positions(res_code)
                except Exception as e:
                    err_msg = "Unable to calculate Kidera: {}. Chain: {}, Residue Number: {}\nSkipping structure due to detected defect".format(
                        e, chain, res_num
                    )
                    logging.warning(
                        "Skipping item %s due to detected defect", self.pdb_code
                    )
                    skip_item = True
                    # not storing under defects any more
                    # self.issues["defects"][self.pdb_code] = {
                    #     "reason": err_msg,
                    #     "physio": self.is_physio,
                    # }
                    break

                kid_values_list = kid_dframe.loc[res_code].to_list()
                res_data["kidera"] = kid_values_list

                # add partial sum of kidera factors
                kidera_sum = [x + y for x, y in zip(kidera_sum, res_data["kidera"])]

                # print(res_data)
                # use pdbe.encode_range instead to account for gaps in 3D model comparing to the sequence (see 3pyi_1.pdb)
                # V2 only: encode_range = (max(1,res_num-3), min(res_num+3, len_seq))  # holds the encode range for the residue r as (r-3, r+3) within the boundaries of the sequence (starting with 1)
                encode_range = (
                    max(1, res_num - 3),
                    res_num + 3,
                )  # holds the encode range for the residue r as (r-3, r+3) starting with 1
                res_data["encode_range"] = encode_range
                encode_range_list = list(range(encode_range[0], encode_range[1] + 1, 1))
                # from_full_seq = encoded_chain[encode_range[0]-1:encode_range[1]-1]
                try:
                    from_range = self.pdbe.encode_range(
                        chain=chain, range=encode_range_list
                    )
                    res_data["pdbencode"] = from_range
                except Exception as e:
                    err_msg = "Unable to PDBEncode: {}. Chain: {}, Residue Number: {}, Encode Range: {}\nSkipping structure due to detected defect".format(
                        e, chain, res_num, encode_range_list
                    )
                    logging.warning(
                        "Skipping item %s due to detected defect", self.pdb_code
                    )
                    skip_item = True
                    # not storing under defects any more
                    # self.issues["defects"][self.pdb_code] = {
                    #     "reason": err_msg,
                    #     "physio": self.is_physio,
                    # }
                    break

                # count residue types
                unknown = True
                if res_data["res_name"] in res_hydrophobic:
                    chain_res_type_count["num_res_hydrophobic"] += 1
                    unknown = False
                if res_data["res_name"] in res_hydrophilic:
                    chain_res_type_count["num_res_hydrophilic"] += 1
                    unknown = False
                if res_data["res_name"] in res_negative:
                    chain_res_type_count["num_res_negative"] += 1
                    unknown = False
                if res_data["res_name"] in res_positive:
                    chain_res_type_count["num_res_positive"] += 1
                    unknown = False
                if res_data["res_name"] in res_aromatic:
                    chain_res_type_count["num_res_aromatic"] += 1
                    unknown = False
                if unknown:
                    chain_res_type_count["num_res_unknown"] += 1
                chain_pdbencode += res_data["pdbencode"]

                if store_residue_data:
                    chains[chain]["residues"][res_num] = res_data

            if skip_item:
                break
            else:
                # calculate average sum of kidera factors
                kidera_avg = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
                if num_chain_residues > 0:
                    kidera_avg = [x / num_chain_residues for x in kidera_sum]
                chains[chain]["kidera"] = kidera_avg
                chains[chain]["num_chain_residues"] = num_chain_residues
                chains[chain]["chain_res_type_count"] = chain_res_type_count
                chains[chain]["pdbencode"] = chain_pdbencode
                # aggregate chain values into interface
                num_group_residues += chains[chain]["num_chain_residues"]
                group_res_type_count["num_res_hydrophobic"] += chains[chain][
                    "chain_res_type_count"
                ]["num_res_hydrophobic"]
                group_res_type_count["num_res_hydrophilic"] += chains[chain][
                    "chain_res_type_count"
                ]["num_res_hydrophilic"]
                group_res_type_count["num_res_negative"] += chains[chain][
                    "chain_res_type_count"
                ]["num_res_negative"]
                group_res_type_count["num_res_positive"] += chains[chain][
                    "chain_res_type_count"
                ]["num_res_positive"]
                group_res_type_count["num_res_aromatic"] += chains[chain][
                    "chain_res_type_count"
                ]["num_res_aromatic"]
                group_res_type_count["num_res_unknown"] += chains[chain][
                    "chain_res_type_count"
                ]["num_res_unknown"]
                group_pdbencode += chain_pdbencode

            # interface_residues[chain]["pdbencode"] = encoded_chain
        if skip_item:
            return dict(), True, err_msg

        # calculate frequencies of residue types at the interface
        group_res_type_freq = {
            "hydrophobic": 0,
            "hydrophilic": 0,
            "negative": 0,
            "positive": 0,
            "aromatic": 0,
            "unknown": 0,
        }
        if num_group_residues > 0:
            group_res_type_freq["hydrophobic"] = (
                group_res_type_count["num_res_hydrophobic"] / num_group_residues
            )
            group_res_type_freq["hydrophilic"] = (
                group_res_type_count["num_res_hydrophilic"] / num_group_residues
            )
            group_res_type_freq["negative"] = (
                group_res_type_count["num_res_negative"] / num_group_residues
            )
            group_res_type_freq["positive"] = (
                group_res_type_count["num_res_positive"] / num_group_residues
            )
            group_res_type_freq["aromatic"] = (
                group_res_type_count["num_res_aromatic"] / num_group_residues
            )
            group_res_type_freq["unknown"] = (
                group_res_type_count["num_res_unknown"] / num_group_residues
            )

        # count chars and frequency of the structural alphabet at the interface
        group_pdbencode_count = self._count_chars_pdbencode(group_pdbencode)
        group_pdbencode_freq = self._freq_chars_pdbencode(group_pdbencode)

        entry = {
            "chains": chains,
            "num_group_residues": num_group_residues,
            "group_res_type_count": group_res_type_count,
            "group_res_type_freq": group_res_type_freq,
            "group_pdbencode": group_pdbencode,
            "group_pdbencode_count": group_pdbencode_count,
            "group_pdbencode_freq": group_pdbencode_freq,
        }

        return entry, False, err_msg

    def _count_chars_pdbencode(self, encoded_sequence: str) -> dict:
        pdbencode_alphabet = list(string.ascii_uppercase)
        pdbencode_alphabet.remove("Z")  # Z does not belong to the PDBEncode alphabet
        count_dict = dict([(key, 0) for key in pdbencode_alphabet])
        for char in encoded_sequence:
            count_dict[char] += 1
        return count_dict

    def _freq_chars_pdbencode(self, encoded_sequence: str) -> dict:
        len_seq = len(encoded_sequence)
        pdbencode_alphabet = list(string.ascii_uppercase)
        pdbencode_alphabet.remove("Z")  # Z does not belong to the PDBEncode alphabet
        count_dict = dict([(key, 0) for key in pdbencode_alphabet])
        freq_dict = dict([(key, 0) for key in pdbencode_alphabet])
        if len_seq > 0:
            for char in encoded_sequence:
                count_dict[char] += 1
            for key in count_dict:
                freq_dict[key] = count_dict[key] / len_seq
        return freq_dict

    def _get_pdb_lab_source(self, pdb_id: str, pdb_file: str) -> str:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", BiopythonWarning)

            parser = PDBParser()
            structure = parser.get_structure(pdb_id, pdb_file)

            header = structure.header["head"]
            if header.endswith("xxxx"):
                # Levy file
                return "levy"
            else:
                # Dunbrack file
                return "dunbrack"

    def _output_chain_residue_to_csv(self, region_data: dict, pdb_id: str):
        res_id_keys = set()
        chain_id_keys = set()

        for chain in self.pops_residues:
            chain_id = pdb_id + "_" + chain

            if chain_id not in chain_id_keys:
                chain_dict = {
                    "chain_id": chain_id,
                    "name": chain,
                    "core_numres": region_data["core"]["chains"][chain][
                        "num_chain_residues"
                    ],
                    "core_k1": region_data["core"]["chains"][chain]["kidera"][0],
                    "core_k2": region_data["core"]["chains"][chain]["kidera"][1],
                    "core_k3": region_data["core"]["chains"][chain]["kidera"][2],
                    "core_k4": region_data["core"]["chains"][chain]["kidera"][3],
                    "sface_numres": region_data["surface"]["chains"][chain][
                        "num_chain_residues"
                    ],
                    "sface_k1": region_data["surface"]["chains"][chain]["kidera"][0],
                    "sface_k2": region_data["surface"]["chains"][chain]["kidera"][1],
                    "sface_k3": region_data["surface"]["chains"][chain]["kidera"][2],
                    "sface_k4": region_data["surface"]["chains"][chain]["kidera"][3],
                    "iface_numres": region_data["interface"]["chains"][chain][
                        "num_chain_residues"
                    ],
                    "iface_k1": region_data["interface"]["chains"][chain]["kidera"][0],
                    "iface_k2": region_data["interface"]["chains"][chain]["kidera"][1],
                    "iface_k3": region_data["interface"]["chains"][chain]["kidera"][2],
                    "iface_k4": region_data["interface"]["chains"][chain]["kidera"][3],
                }
                # CHAIN node
                self.node_chain_csv.writerow(chain_dict.values())
                # PDB->CHAIN relation
                self.rel_pdb_chain_csv.writerow([pdb_id, chain_id, "has"])

                chain_id_keys.add(chain_id)

                for res_label in self.pops_residues[chain]["residues"]:
                    res_id = pdb_id + "_" + chain + "_" + str(res_label)

                    if res_id not in res_id_keys:
                        # RESIDUE node
                        res_data = {
                            "residue-id": res_id,
                            "name": self.pops_residues[chain]["residues"][res_label][
                                "res_name"
                            ],
                            "num": res_label,
                            "q_sasa": self.pops_residues[chain]["residues"][res_label][
                                "q_sasa"
                            ],
                            "location": self.pops_residues[chain]["residues"][
                                res_label
                            ]["location"],
                        }
                        # res_data.update(pops_residues[chain]["residues"][res_label])
                        self.node_residue_csv.writerow(res_data.values())
                        # PDB->RESIDUE relation
                        self.rel_pdb_residue_csv.writerow([pdb_id, res_id, "has"])
                        # CHAIN->RESIDUE relation
                        self.rel_chain_residue_csv.writerow(
                            [chain_id, res_id, "chain_res"]
                        )
                        res_id_keys.add(res_id)
                    else:
                        raise Exception("unexpected duplicate residue key: " + res_id)
            else:
                raise Exception("unexpected duplicate chain key: " + chain_id)

    def _run_get_alphafold_pLDDT(self) -> bool:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", BiopythonWarning)

            parser = PDBParser()
            structure = parser.get_structure(self.pdb_code, self.pdb_file)

            for res_num in self.pops_residues["A"]["residues"]:
                pops_res_name = self.pops_residues["A"]["residues"][res_num]["res_name"]
                structure_res_name = (
                    structure.child_dict[0]
                    .child_dict["A"]
                    .child_list[res_num - 1]
                    .resname
                )
                if pops_res_name != structure_res_name:
                    self.result["processentry_err"] = (
                        "Residue names between pops and structure do not match: {}/{}".format(
                            pops_res_name, structure_res_name
                        )
                    )
                    self.result["processentry_res"] = False
                    return False
                self.pops_residues["A"]["residues"][res_num]["plddt"] = (
                    structure.child_dict[0]
                    .child_dict["A"]
                    .child_list[res_num - 1]
                    .child_dict["CA"]
                    .bfactor
                )

        self.result["processentry_res"] = True
        return True

    def _clean_up(self, pdb_file: str, tmp_pdb_folder: str):
        try:
            if not self.homodimers:
                os.remove(pdb_file)
                for f in os.listdir(tmp_pdb_folder):
                    if f.endswith(".pdb"):
                        os.remove(os.path.join(tmp_pdb_folder, f))
        except Exception:
            pass


def _load_dictionary_from_csv_file(
    file_name: str, file_type: str, dict_key: str, file_not_found_fails=True
) -> dict:
    if not os.path.exists(file_name) and not file_not_found_fails:
        return dict()

    csv_file = ICSV(
        filename=file_name,
        file_type=file_type,
        delimiter=";",
    )
    data_dict = csv_file.get_dict_items(dict_key)
    csv_file.close()
    return data_dict
