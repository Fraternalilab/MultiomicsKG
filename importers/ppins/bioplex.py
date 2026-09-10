from parsers.icsv import ICSV
import importers.ppins.proteinids as pids
import importers.ppins.ppins as ppins


class BioPlexImporter(ICSV):
    """
    Represents a BioPlex importer
    https://bioplex.hms.harvard.edu/interactions.php#datasets
    https://bioplex.hms.harvard.edu/data/BioPlex_interactionList_v4a.tsv

    TODO: Updated BioPlex dataset is available
    """

    def __init__(self, bioplex_filename: str):
        super().__init__(
            filename=bioplex_filename,
            file_type="bioplex",
            max_items=-1,
            delimiter="\t",
            header_filter_prefix="",
            compressed=False
        )

    def importfile(self):
        pass

    def merge_entries(self):
        # performs standard import (calls process_item())
        super().importfile()
        super().close()

    def process_item(self, item: dict):
        protein1_uniacc = item["UniprotA"]
        protein2_uniacc = item["UniprotB"]

        if (
            protein1_uniacc not in pids.uniprot_ids_dict
            or protein2_uniacc not in pids.uniprot_ids_dict
        ):
            # TODO: BioPlex dataset contains many Uniprot accession IDs like Q8N8A2-2 (other isoforms)
            # print(
            #     f"Unexpected interaction bewtween {protein1_uniacc} and {protein2_uniacc}"
            # )
            pass
        else:
            interaction_pair = sorted([protein1_uniacc, protein2_uniacc])
            key = f"{interaction_pair[0]}_{interaction_pair[1]}"

            entry = dict()
            if key in ppins.ppins_dict:
                entry = ppins.ppins_dict[key]
            else:
                entry = ppins.new_entry()
            entry["interaction_pair"] = interaction_pair
            entry["bioplex_p_int"] = float(item["p(Interaction)"])
            entry["bioplex_p_wrong"] = float(item["p(Wrong)"])
            entry["experiment_type"].append("ms")
            entry["sources"].add("bioplex")
            ppins.ppins_dict[key] = entry
