from parsers.icsv import ICSV
import importers.ppins.proteinids as pids
import importers.ppins.ppins as ppins


class CellY2HHIUnionImporter(ICSV):
    """
    Represents a CellY2HHIUnion (vidal20) importer ====> THIS IS THE SAME AS HURI
    https://doi.org/10.1038/s41586-020-2188-x
    http://www.interactome-atlas.org/data/HI-union.tsv

    Updated from an earlier vidal14 dataset:
    https://www.ncbi.nlm.nih.gov/pubmed/25416956
    http://www.interactome-atlas.org/data/HI-II-14.tsv
    """

    def __init__(self, cell_y2h_hi_union_filename: str):
        super().__init__(
            filename=cell_y2h_hi_union_filename,
            file_type="cell_y2h_hi_union",
            max_items=-1,
            delimiter="\t",
            header_filter_prefix="",
            compressed=False,
            fieldnames=["Interactor1", "Interactor2"],
        )

    def importfile(self):
        pass

    def merge_entries(self):
        # performs standard import (calls process_item())
        super().importfile()
        super().close()

    def process_item(self, item: dict):

        if len(item) == 2:
            protein1_uniacc = pids.ensembl_gen_dict.get(item["Interactor1"], None)
            protein2_uniacc = pids.ensembl_gen_dict.get(item["Interactor2"], None)

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
                        # TODO: even within same database there are duplicate interaction pairs
                        entry = ppins.ppins_dict[key]
                    else:
                        entry = ppins.new_entry()
                    entry["interaction_pair"] = interaction_pair
                    entry["experiment_type"].append("y2h")
                    entry["sources"].add("cell_y2h_hi")
                    ppins.ppins_dict[key] = entry
