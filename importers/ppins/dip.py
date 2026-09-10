from parsers.icsv import ICSV
import importers.ppins.proteinids as pids
import importers.ppins.ppins as ppins


class DIPImporter(ICSV):
    """
    Represents a DIP importer
    https://dip.doe-mbi.ucla.edu/dip/Download.cgi (Requires registration)
    Downloaded H.Sapiens dataset Hsapi20170205 from https://dip.doe-mbi.ucla.edu/dip/Download.cgi?SM=7&TX=9606
    https://dip.doe-mbi.ucla.edu/dip/File.cgi?FN=2017/tab25/Hsapi20170205.txt
    """

    def __init__(self, dip_filename: str):
        super().__init__(
            filename=dip_filename,
            file_type="dip",
            max_items=-1,
            delimiter="\t",
            header_filter_prefix="",
            compressed=False,
        )

    def importfile(self):
        pass

    def merge_entries(self):
        self.entries = dict()
        # performs standard import (calls process_item())
        super().importfile()
        super().close()

    def process_item(self, item: dict):
        interactor1 = item["ID interactor A"].split("|uniprotkb:")
        interactor2 = item["ID interactor B"].split("|uniprotkb:")

        if len(interactor1) >= 2 and len(interactor2) >= 2:
            protein1_uniacc = interactor1[1]
            protein2_uniacc = interactor2[1]
            if (
                protein1_uniacc not in pids.uniprot_ids_dict
                or protein2_uniacc not in pids.uniprot_ids_dict
            ):
                # print(
                #     f"Unexpected interaction bewtween {protein1_uniacc} and {protein2_uniacc}, usually due to non-human interactor"
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
                entry["dip_detection_method"] = item[
                    "Interaction detection method(s)"
                ].split("|")
                entry["sources"].add("dip")
                ppins.ppins_dict[key] = entry
        # else:
        #     # unexpected format usually due to missing uniprotkb id
        #     print(f"Unexpected format of interactor in item: {item}")
