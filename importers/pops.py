import parsers.csvfile
import time
import datetime
import csv
import procedures.pops
import logging
from Bio.PDB import PDBList


class POPSImporter:
    """Represents a POPS importer"""

    def __init__(self, filename: str, pdb_database_root: str, scrap_folder: str):
        self.data = []
        self.scrap_folder = scrap_folder
        self.pdb_database_root = pdb_database_root
        self.approx_num_items = 173000
        self.parser = parsers.csvfile.CSVParser(
            filename, "PDB-POPS", self.approx_num_items
        )
        self.num_sel = 0
        self.num_read = 0
        self.num_imported = 0
        self.pops = procedures.pops.POPS(scrap_folder=self.scrap_folder)
        self.report_line = (
            "POPS items imported: {num_items:7d}   |   "
            + "Items/s: {speed:8.2f}   |   "
            + "Time elapsed: {time_elapsed}   |   "
            + "ETA: {eta}"
        )

    def close(self):
        self.parser.close()

    def importdb(self):
        start_counter = time.perf_counter()
        stop_counter = start_counter
        start = time.time()

        node_residue_file = open("node_residue.csv", "w")
        node_residue_csv = csv.writer(node_residue_file, delimiter=";", quotechar="'")
        node_residue_csv.writerow(
            [
                "id:ID(residue-id)",
                "res_name",
                "res_num",
                "icode",
                "q_sasa",
                "phob",
                "phil",
                "surf",
                "loc",
            ]
        )

        rel_pdb_residue_file = open("rel_pdb_residue.csv", "w")
        rel_pdb_residue_csv = csv.writer(
            rel_pdb_residue_file, delimiter=";", quotechar="'"
        )
        rel_pdb_residue_csv.writerow(
            [":START_ID(pdb-id)", ":END_ID(residue-id)", ":TYPE"]
        )

        pdbl = PDBList()

        while True:
            item = self.parser.get_next_item()

            self.num_read += 1

            res_id_keys = set()

            if item is None:
                print("\nNo more items")
                node_residue_file.close()
                rel_pdb_residue_file.close()
                break

            pdb_id = item["id:ID(pdb-id)"]

            logging.info("Processing " + pdb_id)

            # pdbl.retrieve_pdb_file(
            #     pdb_code=pdb_id, file_format="pdb", pdir=self.scrap_folder
            # )

            folder_name = pdb_id[1:3:1]

            pdb_file = (
                self.pdb_database_root + "/" + folder_name + "/pdb" + pdb_id + ".ent"
            )

            pops_residues, msg = self.pops.run(pdb_file, pdb_id)

            if pops_residues is None:
                logging.error(
                    "Unable to POPS the PDB file " + pdb_id + " (skipping): " + msg
                )
                continue

            for chain in pops_residues:
                for res_label in pops_residues[chain]["residues"]:
                    res_id = pdb_id + "_" + chain + "_" + str(res_label)

                    if res_id not in res_id_keys:
                        # RESIDUE node
                        res_data = {"residue-id": res_id}
                        res_data.update(pops_residues[chain]["residues"][res_label])
                        node_residue_csv.writerow(res_data.values())
                        # PDB->RESIDUE relation
                        rel_pdb_residue_csv.writerow([pdb_id, res_id, "has"])

                        res_id_keys.add(res_id)
                    else:
                        raise Exception("unexpected duplicate residue key: " + res_id)

            self.num_imported += 1
            if self.num_read % 1000 == 0:
                stop_counter = time.perf_counter()
                speed = 1000 / (stop_counter - start_counter)
                eta = str(
                    datetime.timedelta(seconds=round(self.approx_num_items / speed))
                )
                time_elapsed = str(
                    datetime.timedelta(seconds=round(time.time() - start))
                )
                print(
                    self.report_line.format(
                        num_items=self.num_imported,
                        speed=speed,
                        time_elapsed=time_elapsed,
                        eta=eta,
                    ),
                    end="\r",
                )
                start_counter = time.perf_counter()
