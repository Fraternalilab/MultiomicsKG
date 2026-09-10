from parsers.icsv import ICSV
import logging
import csv


class SNVResidues(ICSV):
    def __init__(
        self,
        node_snv_filename: str,
        rel_protein_snv_filename: str,
        rel_protein_pdb_filename: str,
        rel_protein_model_filename: str,
        rel_snv_residue_filename: str,
        homodimers=False,
    ):
        super().__init__(
            filename=node_snv_filename,
            file_type="node_snv",
            max_items=-1,
            delimiter=";",
        )

        rel_protein_snv_csv = ICSV(
            rel_protein_snv_filename, "rel_protein_snv", max_items=-1, delimiter=";"
        )
        self.rel_protein_snv_db = rel_protein_snv_csv.get_dict_items(":END_ID(snv-id)")
        rel_protein_snv_csv.close()

        rel_protein_pdb_csv = ICSV(
            rel_protein_pdb_filename, "rel_protein_pdb", max_items=-1, delimiter=";"
        )
        self.rel_protein_pdb_db = rel_protein_pdb_csv.get_dict_items(
            ":START_ID(protein-id)"
        )
        rel_protein_pdb_csv.close()

        rel_protein_model_csv = ICSV(
            rel_protein_model_filename, "rel_protein_model", max_items=-1, delimiter=";"
        )
        self.rel_protein_model_db = rel_protein_model_csv.get_dict_items(
            ":START_ID(protein-id)"
        )
        rel_protein_model_csv.close()
        self.rel_snv_residue_filename = rel_snv_residue_filename
        self.num_snvs = 0
        self.unmapped_snvs_to_pdb = 0
        self.unmapped_snvs_to_model = 0
        self.insertion_code_pdbs = set()
        self.homodimers = homodimers
        logging.getLogger().setLevel(logging.INFO)

    def importfile(self):
        self._set_output_files()

        # performs standard import (calls process_item())
        super().importfile()

        # print PDBs with insertion code
        logging.info(
            "Total number of PDBs with insertion codes: {pdbs}".format(
                pdbs=len(self.insertion_code_pdbs)
            )
        )

    def _set_output_files(self):
        self.rel_snv_residue_file = open(self.rel_snv_residue_filename, "w")
        self.rel_snv_residue_csv = csv.writer(
            self.rel_snv_residue_file, delimiter=";", quotechar="'"
        )
        self.rel_snv_residue_csv.writerow(
            [":START_ID(snv-id)", ":END_ID(residue-id)", ":TYPE"]
        )

    def close(self):
        print(
            "\nTotal SNVs: {num_snvs}, Unmapped to PDB: {unmapped_snvs_to_pdb}, Unmapped to model: {unmapped_snvs_to_model}".format(
                num_snvs=self.num_snvs,
                unmapped_snvs_to_pdb=self.unmapped_snvs_to_pdb,
                unmapped_snvs_to_model=self.unmapped_snvs_to_model,
            )
        )
        self.rel_snv_residue_file.close()
        super().close()

    def process_item(self, item: dict):
        self.num_snvs += 1
        snv_id = item["id:ID(snv-id)"]
        # if item["source"] == "gnomad":
        #     return
        # if snv_id != "787743":
        #     return
        snv_protein_pos = int(item["protein_position:long"])
        # snv_aminoacids = item["amino_acids"]
        for protein_entry in self.rel_protein_snv_db[snv_id]:
            uniprot_id = protein_entry[":START_ID(protein-id)"]
            if uniprot_id not in self.rel_protein_pdb_db:
                self.unmapped_snvs_to_pdb += 1
            else:
                mappings = self.rel_protein_pdb_db[uniprot_id]
                self._output_snv_residue_map(
                    mappings=mappings,
                    snv_protein_pos=snv_protein_pos,
                    snv_id=snv_id,
                    type="pdb",
                )

            if uniprot_id not in self.rel_protein_model_db:
                self.unmapped_snvs_to_model += 1
            else:
                mappings = self.rel_protein_model_db[uniprot_id]
                self._output_snv_residue_map(
                    mappings=mappings,
                    snv_protein_pos=snv_protein_pos,
                    snv_id=snv_id,
                    type="model",
                )

    def _output_snv_residue_map(
        self, mappings: list, snv_protein_pos: int, snv_id: str, type: str
    ):
        for protein_pdb_rel_entry in mappings:
            pdb_id = protein_pdb_rel_entry[":END_ID(" + type + "-id)"]
            chain = protein_pdb_rel_entry["chain"]
            # sp_beg / sp_end = range of positions from uniprot
            # pdb_beg / pdb_end = range of positions from PDB (Author position)
            # res_beg / res_end = range of positions from SEQRES
            try:
                sp_beg = int(protein_pdb_rel_entry["sp_beg"])
                # TODO: improve with range check
                sp_end = int(protein_pdb_rel_entry["sp_end"])
                if snv_protein_pos >= sp_beg and snv_protein_pos <= sp_end:
                    snv_pdb_pos = 0
                    # pdb files in homodimer dataset use SEQRES (res_beg) as residue numbering
                    snv_pdb_pos_rel_res = (snv_protein_pos - sp_beg) + int(
                        protein_pdb_rel_entry["res_beg"]
                    )
                    snv_pdb_pos_rel_pdb = (snv_protein_pos - sp_beg) + int(
                        protein_pdb_rel_entry["pdb_beg"]
                    )

                    if self.homodimers:
                        snv_pdb_pos = snv_pdb_pos_rel_res
                    else:
                        snv_pdb_pos = snv_pdb_pos_rel_pdb

                    # TODO: some sifts mappings for the same chain render the same uniprot_pos -> pdb_pos
                    #       see Q9UI95
                    # print(
                    #     "Aminoacids: {aminoacids}, Protein: {ppos}, PDB: {pdbpos}".format(
                    #         aminoacids=snv_aminoacids,
                    #         ppos=snv_protein_pos,
                    #         pdbpos=snv_pdb_pos,
                    #     )
                    # )

                    res_id = "{pdb_id}_{chain}_{snv_pdb_pos}".format(
                        pdb_id=pdb_id, chain=chain, snv_pdb_pos=snv_pdb_pos
                    )
                    self.rel_snv_residue_csv.writerow([snv_id, res_id, "maps_res"])
                    # print(protein_pdb_rel_entry)
            except ValueError:
                # logging.error("Unable to handle PDB mappings with insertion code for PDB {pdb_id}, Uniprot {uniprot_id}: {error}".format(pdb_id=pdb_id, uniprot_id=uniprot_id, error=ve))
                self.insertion_code_pdbs.add(pdb_id)
