from parsers.icsv import ICSV
import importers.ppins.proteinids as pids
import importers.ppins.ppins as ppins


class STRINGImporter(ICSV):
    """
    Represents a STRING importer
    https://stringdb-downloads.org/download/protein.links.detailed.v12.0/9606.protein.links.detailed.v12.0.txt.gz
    """

    def __init__(self, string_filename: str):
        super().__init__(
            filename=string_filename,
            file_type="string",
            max_items=-1,
            delimiter=" ",
            header_filter_prefix="",
            compressed=False,
        )

    def importfile(self):
        pass

    def merge_entries(self):
        # performs standard import (calls process_item())
        super().importfile()
        super().close()

    def process_item(self, item: dict):

        protein1 = item["protein1"].split(".")[1]
        protein2 = item["protein2"].split(".")[1]
        protein1_uniacc = pids.ensembl_pro_dict.get(protein1, None)
        protein2_uniacc = pids.ensembl_pro_dict.get(protein2, None)
        if (protein1_uniacc is not None) and (protein2_uniacc is not None):
            if (
                protein1_uniacc not in pids.uniprot_ids_dict
                or protein2_uniacc not in pids.uniprot_ids_dict
            ):
                print(
                    f"Unexpected interaction bewtween {protein1_uniacc} and {protein2_uniacc}"
                )
            else:
                interaction_pair = sorted([protein1_uniacc, protein2_uniacc])
                key = f"{interaction_pair[0]}_{interaction_pair[1]}"

                entry = dict()
                if key in ppins.ppins_dict:
                    # TODO: even within STRING database there are duplicate interaction pairs with different scores
                    entry = ppins.ppins_dict[key]
                else:
                    entry = ppins.new_entry()
                entry["interaction_pair"] = interaction_pair
                # TODO: Filter here only experimental score > 700 (confidence > 0.7, "high-confidence")?
                entry["string_experimental_score"] = float(item["experimental"])
                entry["string_combined_score"] = float(item["combined_score"])
                entry["sources"].add("string")
                ppins.ppins_dict[key] = entry
        # else:
        #     print(f"Item missing at least one interactor's Uniprot Accession: {item}")
