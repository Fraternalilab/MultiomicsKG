from parsers.icsv import ICSV
import importers.ppins.proteinids as pids
import importers.ppins.ppins as ppins


class BioGRIDImporter(ICSV):
    """
    Represents a BioGRID importer
    https://downloads.thebiogrid.org/BioGRID/Release-Archive
    https://downloads.thebiogrid.org/File/BioGRID/Release-Archive/BIOGRID-4.4.236/BIOGRID-ALL-4.4.236.tab3.zip
    """

    def __init__(self, biogrid_filename: str):
        super().__init__(
            filename=biogrid_filename,
            file_type="biogrid",
            max_items=-1,
            delimiter="\t",
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

        protein1_uniacc = item["SWISS-PROT Accessions Interactor A"]
        protein2_uniacc = item["SWISS-PROT Accessions Interactor B"]

        # validate IDs (skip "-" or ids containing "|" )
        if (
            protein1_uniacc == "-"
            or protein2_uniacc == "-"
            or "|" in protein1_uniacc
            or "|" in protein2_uniacc
        ):
            # TODO: Handle cases where [Q9H7B7 - P01892|P30443|P30450|P04439] for example
            # print(
            #     f"Skipping intractable iteration [{protein1_uniacc} - {protein2_uniacc}]"
            # )
            return
        else:
            if (
                protein1_uniacc not in pids.uniprot_ids_dict
                or protein2_uniacc not in pids.uniprot_ids_dict
            ):
                # TODO: Cases found in BioGrid dataset
                # print(
                #     f"Unexpected interaction bewtween {protein1_uniacc} and {protein2_uniacc}"
                # )
                pass
            else:
                interaction_pair = sorted([protein1_uniacc, protein2_uniacc])
                key = f"{interaction_pair[0]}_{interaction_pair[1]}"

                entry = dict()
                if key in ppins.ppins_dict:
                    # TODO: even within BioGrid database there are duplicate interaction pairs
                    entry = ppins.ppins_dict[key]
                else:
                    entry = ppins.new_entry()
                entry["interaction_pair"] = interaction_pair
                entry["biogrid_experimental_system"].append(
                    str.lower(item["Experimental System"])
                )
                entry["biogrid_experimental_system_type"].append(
                    str.lower(item["Experimental System Type"])
                )
                entry["biogrid_throughput"].append(str.lower(item["Throughput"]))
                score = item["Score"]
                if score != "-":
                    entry["biogrid_score"].append(float(score))
                entry["sources"].add("biogrid")
                publication_source = item["Publication Source"]
                entry["publication"].append(publication_source)
                match publication_source:
                    case "PUBMED:22939629":
                        entry["experiment_type"].append("ms")
                        entry["sources"].add("bpms")
                    case "PUBMED:26496610":
                        entry["experiment_type"].append("ms")
                        entry["sources"].add("hein_ms")
                    case "PUBMED:32296183":
                        entry["experiment_type"].append("y2h")
                        entry["sources"].add("huri")
                ppins.ppins_dict[key] = entry
