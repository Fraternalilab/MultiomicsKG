from parsers.icsv import ICSV
import importers.ppins.proteinids as pids
import importers.ppins.ppins as ppins


class MetazoanImporter(ICSV):
    """
    Represents a Metazoan importer

    http://metazoa.med.utoronto.ca/projected_ppi/Hsapiens_projected_PPI_metazoan_conserved_map.txt

    """

    def __init__(self, metazoan_filename: str):
        super().__init__(
            filename=metazoan_filename,
            file_type="metazoan",
            max_items=-1,
            delimiter="\t",
            compressed=False,
            # must skip header and redefine it because of duplicate columns in original source file
            skip_rows=1,
            fieldnames=[
                "Hsapiens Ensembl Gene1",
                "Hsapiens Ensembl Gene2",
                "Hsapiens Gene1",
                "Hsapiens Gene2",
                "Score",
                "Metazoan ComplexID1",
                "Metazoan ComplexID2",
                "Description1",
                "Description2",
            ],
        )

    def importfile(self):
        pass

    def merge_entries(self):
        # performs standard import (calls process_item())
        super().importfile()
        super().close()

    def process_item(self, item: dict):

        if len(item) > 2:
            protein1 = item["Hsapiens Ensembl Gene1"].split(
                "||"
            )  # Hsapiens Gene1 (duplicate column)
            protein2 = item["Hsapiens Ensembl Gene2"].split(
                "||"
            )  # Hsapiens Gene2 (duplicate column)

            for iso_protein1 in protein1:
                iso_protein1_uniacc = pids.ensembl_gen_dict.get(iso_protein1, None)
                if iso_protein1_uniacc is not None:
                    for iso_protein2 in protein2:
                        iso_protein2_uniacc = pids.ensembl_gen_dict.get(
                            iso_protein2, None
                        )
                        if iso_protein2_uniacc is not None:
                            if (
                                iso_protein1_uniacc not in pids.uniprot_ids_dict
                                or iso_protein2_uniacc not in pids.uniprot_ids_dict
                            ):
                                print(
                                    f"Unexpected interaction bewtween {iso_protein1_uniacc} and {iso_protein2_uniacc}"
                                )
                            else:
                                interaction_pair = sorted(
                                    [iso_protein1_uniacc, iso_protein2_uniacc]
                                )
                                key = f"{interaction_pair[0]}_{interaction_pair[1]}"

                                entry = dict()
                                if key in ppins.ppins_dict:
                                    entry = ppins.ppins_dict[key]
                                else:
                                    entry = ppins.new_entry()
                                entry["interaction_pair"] = interaction_pair
                                # TODO: Filter here only experimental score > 700 (confidence > 0.7, "high-confidence")?
                                entry["metazoan_score"] = float(item["Score"])
                                entry["experiment_type"].append("ms")
                                entry["sources"].add("metazoan")
                                ppins.ppins_dict[key] = entry
