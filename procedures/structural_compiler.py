from parsers.icsv import ICSV
import logging
import json
import os
import procedures.common as common
import pickle as pickle
import math

data_results_filename_suffix = ".data_results.pickle"

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


class StructuralCompilerRunner:
    def __init__(
        self,
        node_pdb_filename: str,
        node_model_filename: str,
        scrap_folder: str,
        pdb_issues_filename: str,
        model_issues_filename: str,
        node_chain_filename: str,
        node_residue_filename: str,
        rel_pdb_chain_filename: str,
        rel_model_chain_filename: str,
        rel_pdb_residue_filename: str,
        rel_model_residue_filename: str,
        rel_chain_residue_filename: str,
        rel_residue_residue_filename: str,
        skip_defects: bool,
    ):
        self.nodes_filenames = {
            "pdb": node_pdb_filename,
            "model": node_model_filename,
        }
        self.scrap_folder = scrap_folder
        self.issues_filenames = {
            "pdb": pdb_issues_filename,
            "model": model_issues_filename,
        }
        self.items_read = 0
        self.skip_defects = skip_defects

        self.node_chain_filename = node_chain_filename
        self.node_residue_filename = node_residue_filename
        self.rel_to_chain_filenames = {
            "pdb": rel_pdb_chain_filename,
            "model": rel_model_chain_filename,
        }
        self.rel_to_residue_filenames = {
            "pdb": rel_pdb_residue_filename,
            "model": rel_model_residue_filename,
        }
        self.rel_chain_residue_filename = rel_chain_residue_filename
        self.rel_residue_residue_filename = rel_residue_residue_filename

        self.issues = dict()
        self.result = _make_task_result()

    def run(self):
        self._load_structure_nodes(types=["pdb", "model"])

        self._set_output_files()

        self._load_issues_from_file(types=["pdb", "model"])

        self._select_dispatch_items(types=["pdb", "model"])

        self._close()

    def _load_structure_nodes(self, types: list):
        self.nodes = dict()
        for type in types:
            self.nodes[type] = _load_dictionary_from_csv_file(
                file_name=self.nodes_filenames[type],
                file_type="node_" + type,
                dict_key="id:ID(" + type + "-id)",
            )

    def _load_issues_from_file(self, types: list):
        self.issues = dict()
        # self.issues["items_read"] = 0
        # self.issues["max_items"] = -1
        # self.issues["report"] = dict()
        # read issues files
        for type in types:
            try:
                issues_file = open(self.issues_filenames[type], mode="r")
                self.issues[type] = json.load(issues_file)
                issues_file.close()
            except json.JSONDecodeError as j:
                logging.fatal(
                    "Error whilst parsing existing issue file '%s': %s",
                    self.issues_filenames[type],
                    j,
                )
                exit(code=1)
            except Exception:
                logging.fatal("Issues file not found: %s", self.issues_filenames[type])
                exit(code=1)

    def _select_dispatch_items(self, types: list):
        for type in types:
            for structure_id in self.nodes[type]:
                self.items_read += 1

                if self.skip_defects and structure_id in self.issues[type]["report"]:
                    if (
                        self.issues[type]["report"][structure_id]["general_res"]
                        is False
                    ):
                        logging.warning(
                            "Skipping entry with known issue: %s", structure_id
                        )
                        continue

                logging.info("Processing %s %s", type, structure_id)

                # scrap_pdb_folder = self.scrap_folder + "/" + structure_id
                subfolder_divider_code = structure_id[1:3]
                scrap_pdb_folder = "{0}/{1}/{2}".format(
                    self.scrap_folder, subfolder_divider_code, structure_id
                )
                data_results_filename = "{0}/{1}{2}".format(
                    scrap_pdb_folder, structure_id, data_results_filename_suffix
                )

                (
                    result,
                    pops_residues,
                    region_data,
                    contacts_list,
                ) = self._load_data_results_file(
                    data_results_filename=data_results_filename
                )

                if result["processentry_res"]:
                    self._output_chain_residue_to_csv(
                        region_data=region_data,
                        pops_residues=pops_residues,
                        structure_id=structure_id,
                        type=type,
                    )
                else:
                    logging.error("Data unavailable for " + structure_id)

                if result["contacts_list_res"]:
                    self._output_residue_residue_to_csv(
                        contacts_list=contacts_list, structure_id=structure_id
                    )
                else:
                    logging.error("Contacts list unavailable for " + structure_id)

    def _set_output_files(self):
        # (
        #     self.pdb_residues_progress_file,
        #     self.pdb_residues_progress,
        # ) = common.create_or_append_file(
        #     self.pdb_residues_progress_filename, ["pdb-id"]
        # )
        self.node_chain_file, self.node_chain_csv = common.create_or_append_file(
            self.node_chain_filename,
            [
                "id:ID(chain-id)",
                "name:string",
                "core_numres:long",
                "core_k1:double",
                "core_k2:double",
                "core_k3:double",
                "core_k4:double",
                "sface_numres:long",
                "sface_k1:double",
                "sface_k2:double",
                "sface_k3:double",
                "sface_k4:double",
                "iface_numres:long",
                "iface_k1:double",
                "iface_k2:double",
                "iface_k3:double",
                "iface_k4:double",
            ],
            force_create=True,
        )
        self.node_residue_file, self.node_residue_csv = common.create_or_append_file(
            self.node_residue_filename,
            [
                "id:ID(residue-id)",
                "name:string",
                "number:long",
                "q_sasa:double",
                "location:string",
                "plddt:double",
            ],
            force_create=True,
        )
        self.rel_to_chain_files = dict()
        self.rel_to_chain_csvs = dict()
        self.rel_to_residue_files = dict()
        self.rel_to_residue_csvs = dict()
        for type in ["pdb", "model"]:
            (
                self.rel_to_chain_files[type],
                self.rel_to_chain_csvs[type],
            ) = common.create_or_append_file(
                self.rel_to_chain_filenames[type],
                [
                    ":START_ID(" + type + "-id)",
                    ":END_ID(chain-id)",
                    ":TYPE",
                ],
                force_create=True,
            )
            (
                self.rel_to_residue_files[type],
                self.rel_to_residue_csvs[type],
            ) = common.create_or_append_file(
                self.rel_to_residue_filenames[type],
                [
                    ":START_ID(" + type + "-id)",
                    ":END_ID(residue-id)",
                    ":TYPE",
                ],
                force_create=True,
            )
        (
            self.rel_chain_residue_file,
            self.rel_chain_residue_csv,
        ) = common.create_or_append_file(
            self.rel_chain_residue_filename,
            [
                ":START_ID(chain-id)",
                ":END_ID(residue-id)",
                ":TYPE",
            ],
            force_create=True,
        )

        (
            self.rel_residue_residue_file,
            self.rel_residue_residue_csv,
        ) = common.create_or_append_file(
            self.rel_residue_residue_filename,
            [
                ":START_ID(residue-id)",
                ":END_ID(residue-id)",
                ":TYPE",
            ],
            force_create=True,
        )

    def _close(self):
        # self.pdb_residues_progress_file.close()
        self.node_chain_file.close()
        self.node_residue_file.close()
        for type in ["pdb", "model"]:
            self.rel_to_chain_files[type].close()
            self.rel_to_residue_files[type].close()
        self.rel_chain_residue_file.close()
        self.rel_residue_residue_file.close()

    def _load_data_results_file(
        self, data_results_filename: str
    ) -> tuple[dict, dict, dict, list]:
        pops_residues = dict()
        region_data = dict()
        contacts_list = list()
        result = _make_task_result()
        if os.path.exists(data_results_filename):
            with open(data_results_filename, "rb") as file:
                data_results = pickle.load(file)
                if "results" in data_results:
                    result = data_results["results"]
                if "pops_residues" in data_results:
                    pops_residues = data_results["pops_residues"]
                if "region_data" in data_results:
                    region_data = data_results["region_data"]
                if "contacts_list" in data_results:
                    contacts_list = data_results["contacts_list"]

        return result, pops_residues, region_data, contacts_list

    def _output_chain_residue_to_csv(
        self, region_data: dict, pops_residues: dict, structure_id: str, type: str
    ):
        res_id_keys = set()
        chain_id_keys = set()

        for chain in pops_residues:
            chain_id = structure_id + "_" + chain

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
                # STRUCTURE->CHAIN relation
                self.rel_to_chain_csvs[type].writerow([structure_id, chain_id, "has"])

                chain_id_keys.add(chain_id)

                for res_label in pops_residues[chain]["residues"]:
                    res_id = structure_id + "_" + chain + "_" + str(res_label)

                    if res_id not in res_id_keys:
                        # RESIDUE node
                        res_data = dict()

                        res_data["residue-id"] = res_id
                        res_data["name"] = pops_residues[chain]["residues"][res_label][
                            "res_name"
                        ]
                        res_data["number"] = res_label
                        res_data["q_sasa"] = ""
                        if not math.isnan(
                            pops_residues[chain]["residues"][res_label]["q_sasa"]
                        ):
                            res_data["q_sasa"] = float(
                                pops_residues[chain]["residues"][res_label]["q_sasa"]
                            )
                        res_data["location"] = pops_residues[chain]["residues"][
                            res_label
                        ]["location"]
                        # add plddt if residue comes from AlphaFold model
                        if type == "model":
                            res_data["plddt"] = float(
                                pops_residues[chain]["residues"][res_label]["plddt"]
                            )
                        else:
                            res_data["plddt"] = 100.0
                        # res_data.update(pops_residues[chain]["residues"][res_label])
                        self.node_residue_csv.writerow(res_data.values())
                        # STRUCTURE->RESIDUE relation
                        self.rel_to_residue_csvs[type].writerow(
                            [structure_id, res_id, "has"]
                        )
                        # CHAIN->RESIDUE relation
                        self.rel_chain_residue_csv.writerow(
                            [chain_id, res_id, "chain_res"]
                        )
                        res_id_keys.add(res_id)
                    else:
                        raise Exception("unexpected duplicate residue key: " + res_id)
            else:
                raise Exception("unexpected duplicate chain key: " + chain_id)

    def _output_residue_residue_to_csv(self, contacts_list: list, structure_id: str):
        res_res_pairs = set()

        for contact_pair in contacts_list:
            resA = contact_pair[0]
            resB = contact_pair[1]
            keyA = resA + "-" + resB
            keyB = resB + "-" + resA
            if keyA not in res_res_pairs and keyB not in res_res_pairs:
                self.rel_residue_residue_csv.writerow([resA, resB, "contact"])
                res_res_pairs.add(keyA)
            else:
                logging.warning(
                    f"Duplicate contact pair: {resA} / {resB} ({structure_id})"
                )


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
