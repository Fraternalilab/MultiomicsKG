import logging
from rpy2.robjects.packages import importr
import rpy2.rinterface_lib.callbacks as r_cbs
import subprocess
import os
import shutil
import csv
import xmlrpc.client as xmlrpclib
import argparse

# Redirect console output from R to generic logging
# Ref: https://rpy2.github.io/doc/latest/html/callbacks.html

logger = logging.getLogger()


def r_consolewrite_print_cb(s: str):
    logger.info(s)


def r_consolewrite_warnerror_cb(s: str):
    logger.warning(s)


r_cbs.consolewrite_print = r_consolewrite_print_cb
r_cbs.consolewrite_warnerror = r_consolewrite_warnerror_cb


class POPScomp:
    """Class that interacts with R to run POPScomp"""

    def __init__(self, scrap_folder: str, visualise: bool):
        self.scrap_folder = scrap_folder
        self.popsr = importr("POPSR")
        self.rprobe = "1.4"  # default value for rprobe
        self.delta_file_name = "deltaSASA.rpopsResidue"
        self.iso_file_name = "isoSASA.rpopsResidue"
        self.id_file_name = "id.rpopsResidue"
        self.interface_dsasa_threshold = 15.0
        self.core_q_sasa_threshold = 0.15
        self.visualise = visualise

    def run(self, pdb_file_path: str, pdb_code: str):
        pops_residues = dict()
        err_msg = ""
        try:
            # copy pdb_file to the scrap folder
            pdb_file_name = pdb_code + ".pdb"
            # non-divided folder
            # scrap_folder = self.scrap_folder + "/" + pdb_code
            subfolder_divider_code = pdb_code[1:3]
            scrap_folder = "{0}/{1}/{2}".format(
                self.scrap_folder, subfolder_divider_code, pdb_code
            )
            if not os.path.exists(scrap_folder):
                os.makedirs(scrap_folder)
            tmp_pdb_file_path = scrap_folder + "/" + pdb_file_name
            shutil.copyfile(pdb_file_path, tmp_pdb_file_path)

            # run first pass of pops to generate id files
            args = [
                "pops",
                "--outDirName",
                scrap_folder,
                "--rout",
                "--atomOut",
                "--residueOut",
                "--chainOut",
                "--rProbe",
                self.rprobe,
                "--pdb",
                tmp_pdb_file_path,
                "1> POPScomp.o 2> POPScomp.e",
            ]
            result = subprocess.run(
                args=args, check=True, text=True, capture_output=True
            )
            if result.returncode != 0:
                err_msg = "POPS first-pass raised an error: \n {error}".format(
                    error=result.stderr
                )
                logging.error(err_msg)
                return None, err_msg

            if (
                len(result.stderr) > 0
            ):  # pops does not return error code when error happens (system error = 0)
                err_msg = "POPS first-pass raised a warning: \n {error}".format(
                    error=result.stderr
                )
                logging.error(err_msg)

            # run popscomp to generate delta files
            self.popsr.popscompR(inputPDB=pdb_file_name, outDir=scrap_folder)

            # read iso file (to find final core and surface/interface residues)
            iso_filename_path = scrap_folder + "/" + self.iso_file_name
            id_filename_path = scrap_folder + "/" + self.id_file_name
            baseline_filename_path = ""
            using_id_file = False
            if os.path.isfile(iso_filename_path):
                baseline_filename_path = iso_filename_path
            else:
                using_id_file = True
                logging.info(
                    "ISO file not found. Using ID file, probably single-chain PDB"
                )
                baseline_filename_path = id_filename_path

            with open(file=baseline_filename_path, mode="r") as iso_file:
                field_names = [
                    "ResidNe",
                    "Chain",
                    "ResidNr",
                    "iCode",
                    "Phob/A^2",
                    "Phil/A^2",
                    "SASA/A^2",
                    "Q(SASA)",
                    "N(overl)",
                    "Surf/A^2",
                ]
                iso_reader = csv.DictReader(
                    iso_file, delimiter="\t", fieldnames=field_names
                )
                try:
                    for item in iso_reader:
                        # skip the header line
                        if iso_reader.line_num == 1:
                            continue
                        q_sasa = float(item["Q(SASA)"].replace(" ", ""))  # remove space
                        chain = item["Chain"].replace(" ", "")
                        res_num = int(item["ResidNr"].replace(" ", ""))
                        residue = {
                            "res_name": item["ResidNe"].replace(" ", ""),
                            "q_sasa": q_sasa,
                        }

                        if chain not in pops_residues:
                            # create a dictionary for residues if the chain does not exist
                            pops_residues[chain] = {"residues": dict()}

                        # check if residue has been seen before
                        if res_num in pops_residues[chain]["residues"]:
                            err_msg = "{err_msg}\nPOPScomp warning:\nDuplicate residue seen in ISO file: chain={chain}, residue_number={res_num}".format(
                                err_msg=err_msg, chain=chain, res_num=res_num
                            )
                            return None, err_msg

                        if q_sasa <= self.core_q_sasa_threshold:
                            # core residue
                            residue["location"] = "core"
                        else:
                            # surface residue (or potentially interface to be overwritten if found in delta file)
                            residue["location"] = "surface"
                        pops_residues[chain]["residues"][res_num] = residue
                except csv.Error as e:
                    err_msg = "{err_msg}\nPOPScomp error:\nError parsing iso file, line {line}: {error}".format(
                        err_msg=err_msg, line=iso_reader.line_num, error=e
                    )
                    return None, err_msg

            if using_id_file:
                return pops_residues, err_msg

            # read delta file ()
            with open(
                file=scrap_folder + "/" + self.delta_file_name, mode="r"
            ) as delta_file:
                field_names = [
                    "FileLabel",
                    "ResidNe",
                    "Chain",
                    "ResidNr",
                    "iCode",
                    "D_Phob.A.2",
                    "D_Phil.A.2",
                    "D_SASA.A.2",
                ]
                delta_reader = csv.DictReader(
                    delta_file, delimiter=" ", fieldnames=field_names
                )
                try:
                    for item in delta_reader:
                        # skip the header line
                        if delta_reader.line_num == 1:
                            continue
                        d_sasa_a2 = float(item["D_SASA.A.2"])
                        chain = item["Chain"]
                        res_num = int(item["ResidNr"])

                        # check if chain has been seen in ISO file
                        if chain not in pops_residues:
                            err_msg = "{err_msg}\nPOPScomp warning:\nChain found in delta file not seen in ISO file: chain={chain}, residue_number={res_num}".format(
                                err_msg=err_msg, chain=chain, res_num=res_num
                            )
                            return None, err_msg

                        # residue found in delta file is expected to already exist in ISO
                        if res_num in pops_residues[chain]["residues"]:
                            # update existing entry from ISO file
                            pops_residues[chain]["residues"][res_num][
                                "dsasaa2"
                            ] = d_sasa_a2
                            if d_sasa_a2 > self.interface_dsasa_threshold:
                                pops_residues[chain]["residues"][res_num][
                                    "location"
                                ] = "interface"
                        else:
                            err_msg = "{err_msg}\nPOPScomp warning:\nResidue found in delta file not seen in ISO file: chain={chain}, residue_number={res_num}".format(
                                err_msg=err_msg, chain=chain, res_num=res_num
                            )
                            return None, err_msg
                except csv.Error as e:
                    err_msg = "{err_msg}\nPOPScomp error:\nError parsing delta file, line {line}: {error}".format(
                        err_msg=err_msg, line=delta_reader.line_num, error=e
                    )
                    return None, err_msg
                    # raise Exception('Error parsing delta file, line {}: {}'.format(delta_reader.line_num, e))

            # visualise on pymol
            if self.visualise:
                try:
                    pymolrpc = xmlrpclib.ServerProxy("http://localhost:9123")
                    pymolrpc.delete("all")  # clean up view
                    pymolrpc.load(tmp_pdb_file_path)  # load complex pdb file
                    pymolrpc.hide("cartoon", "all")  # remove default cartoon view
                    pymolrpc.show("ribbon", "all")  # set ribbon view
                    colors_chains = (
                        "magenta",
                        "teal",
                        "blue",
                    )  # list of possible colours to assign to chains (only two really needed)
                    colors_interfaces = ("yellow", "green", "cyan")
                    idx_color = 0
                    for chain in pops_residues:
                        list_res = []  # list of residues on the interface of this chain
                        for res_num in pops_residues[chain]["residues"]:
                            if (
                                pops_residues[chain]["residues"][res_num]["location"]
                                == "interface"
                            ):
                                list_res.append(str(res_num))
                        pymolrpc.select(
                            "sel_int" + chain,
                            "chain " + chain + " and resi " + ",".join(list_res),
                        )  # create a selection of interface residues of this chain
                        pymolrpc.create("interface" + chain, "sel_int" + chain)
                        pymolrpc.delete("sel_int" + chain)
                        pymolrpc.color(colors_chains[idx_color], "chain " + chain)
                        pymolrpc.color(
                            colors_interfaces[idx_color], "interface" + chain
                        )
                        pymolrpc.show("licorice", "interface" + chain)
                        idx_color = idx_color + 1
                    pymolrpc.zoom("all")
                except ConnectionRefusedError as exception:
                    logging.error(
                        "Unable to connect to PyMOL. Run 'pymol -R' in order to visualise the complex"
                    )

            # self.cleanup(pdb_code)
            return pops_residues, err_msg
        except Exception as e:
            # os.remove(tmp_pdb_file_path)
            # self.cleanup(pdb_code)
            err_msg = "{err_msg}\nGeneral error: {exception}".format(
                err_msg=err_msg, exception=e
            )
            return None, err_msg
            # raise exception

    def cleanup(self, pdb_code: str):
        # clean up the scrap folder
        tmp_files = os.listdir(self.scrap_folder)
        for file in tmp_files:
            if (
                file.startswith(pdb_code)
                or file.startswith("deltaSASA")
                or file.startswith("id.rpops")
                or file.startswith("isoSASA")
                or file.startswith("POPScomp")
            ):
                os.remove(self.scrap_folder + "/" + file)

    def load_pops_residues_from_csv(self, filename: str) -> dict:
        pops_residues = dict()
        with open(filename, "r") as file:
            csv_file = csv.DictReader(
                file, delimiter=";", quotechar="'", quoting=csv.QUOTE_NONNUMERIC
            )
            for item in csv_file:
                chain = item["chain"]
                res_num = int(item["res_num"])
                res_data = {
                    "res_name": item["res_name"],
                    "q_sasa": item["q_sasa"],
                    "location": item["location"],
                }
                if chain not in pops_residues:
                    pops_residues[chain] = {
                        "residues": dict(),
                    }
                pops_residues[chain]["residues"][res_num] = res_data

        return pops_residues

    def save_pops_residues_to_csv(self, filename: str, pops_residues: dict):
        with open(filename, "w") as file:
            csv_file = csv.writer(
                file, delimiter=";", quotechar="'", quoting=csv.QUOTE_NONNUMERIC
            )
            csv_file.writerow(["chain", "res_num", "res_name", "q_sasa", "location"])

            for chain in pops_residues.keys():
                for res_num in pops_residues[chain]["residues"].keys():
                    res_name = pops_residues[chain]["residues"][res_num]["res_name"]
                    q_sasa = pops_residues[chain]["residues"][res_num]["q_sasa"]
                    location = pops_residues[chain]["residues"][res_num]["location"]
                    csv_file.writerow([chain, res_num, res_name, q_sasa, location])


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='Run POPScomp on a PDB file and visualise interface residues on PyMOL. Requires PyMOL running in server mode: "pymol -R"'
    )
    parser.add_argument(
        "--pdb_file",
        dest="pdb_file",
        type=str,
        required=True,
        help="path to the PDB file",
    )
    parser.add_argument(
        "--pdb_code", dest="pdb_code", type=str, required=True, help="PDB code"
    )
    parser.add_argument(
        "--scrap",
        dest="temp_folder",
        type=str,
        required=True,
        help="empty folder to store temporary files (all files are removed!!)",
    )

    args = parser.parse_args()

    popscomp = POPScomp(scrap_folder=args.temp_folder, visualise=True)
    popscomp.run(args.pdb_file, args.pdb_code)
