from parsers.icsv import ICSV
import csv


class SIFTSImporter(ICSV):
    """Represents a SIFTS importer"""

    def __init__(
        self,
        filename: str,
        node_protein_filename: str,
        node_pdb_filename: str,
        rel_protein_pdb_filename: str,
        filter_by_node_protein: bool,
        compressed: bool = False,
    ):
        super().__init__(
            filename=filename,
            file_type="sifts",
            max_items=-1,
            delimiter=",",
            header_filter_prefix="#",
            compressed=compressed,
        )
        self.num_read = 0
        self.num_sel = 0
        self.num_imported = 0
        self.node_protein_filename = node_protein_filename
        self.node_pdb_filename = node_pdb_filename
        self.rel_protein_pdb_filename = rel_protein_pdb_filename
        self.protein_set = set()
        self.filter_by_node_protein = filter_by_node_protein

    def close(self):
        print(
            "Total entries: {num_read}, selected: {selected}, mapped: {mapped}".format(
                num_read=self.num_read,
                selected=self.num_sel,
                mapped=self.num_imported,
            )
        )
        self.node_pdb_file.close()
        self.rel_protein_pdb_file.close()

        super().close()

    def _load_proteins_set(self):
        node_protein_csv = ICSV(
            filename=self.node_protein_filename,
            file_type="node_protein",
            max_items=-1,
            delimiter=";",
        )
        # list of chars to be replaced in multi-valued accession_ids
        replace_chars = {"'": "", "[": "", "]": "", " ": ""}
        while True:
            item = node_protein_csv.get_next_item()
            if item is None:
                node_protein_csv.close()
                return

            primary_id = item["accession:ID(protein-id)"]
            accession_ids = item["alternative_accessions:string[]"]
            # remove '[] and space
            accession_ids = "".join(replace_chars.get(s, s) for s in accession_ids)
            accession_ids = accession_ids.split(",")
            self.protein_set.add(primary_id)
            self.protein_set.update(accession_ids)

    def _set_output_files(self):
        self.node_pdb_file = open(self.node_pdb_filename, "w")
        self.node_pdb_csv = csv.writer(self.node_pdb_file, delimiter=";", quotechar="'")
        self.node_pdb_csv.writerow(["id:ID(pdb-id)"])

        self.rel_protein_pdb_file = open(self.rel_protein_pdb_filename, "w")
        self.rel_protein_pdb_csv = csv.writer(
            self.rel_protein_pdb_file, delimiter=";", quotechar="'"
        )
        self.rel_protein_pdb_csv.writerow(
            [
                ":START_ID(protein-id)",
                "chain",
                "res_beg",
                "res_end",
                "pdb_beg",
                "pdb_end",
                "sp_beg",
                "sp_end",
                ":END_ID(pdb-id)",
                ":TYPE",
            ]
        )

    def importfile(self):
        if self.filter_by_node_protein:
            self._load_proteins_set()

        self._set_output_files()

        self.pdb_keys = set()

        # performs standard import (calls process_item())
        super().importfile()

    def process_item(self, item: dict):
        self.num_read += 1

        pdb_id = item["PDB"]
        uniprot_id = item["SP_PRIMARY"]

        if self.filter_by_node_protein and uniprot_id not in self.protein_set:
            return

        self.num_sel += 1
        map_protein_pdb = dict()

        map_protein_pdb["uniprot"] = uniprot_id
        map_protein_pdb["chain"] = item["CHAIN"]
        map_protein_pdb["res_beg"] = item["RES_BEG"]
        map_protein_pdb["res_end"] = item["RES_END"]
        map_protein_pdb["pdb_beg"] = item["PDB_BEG"]
        map_protein_pdb["pdb_end"] = item["PDB_END"]
        map_protein_pdb["sp_beg"] = item["SP_BEG"]
        map_protein_pdb["sp_end"] = item["SP_END"]
        map_protein_pdb["pdb"] = pdb_id
        map_protein_pdb["type"] = "maps"

        # PDB node
        if pdb_id not in self.pdb_keys:
            self.node_pdb_csv.writerow([pdb_id])
            self.pdb_keys.add(pdb_id)

        # PROTEIN->PDB relation
        self.rel_protein_pdb_csv.writerow(map_protein_pdb.values())

        self.num_imported += 1
