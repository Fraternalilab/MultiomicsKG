from parsers.icsv import ICSV
from importers.ppins.biogrid import BioGRIDImporter
from importers.ppins.bioplex import BioPlexImporter
from importers.ppins.cell_y2h_hi_union import CellY2HHIUnionImporter
from importers.ppins.dip import DIPImporter
from importers.ppins.intact import IntActImporter
from importers.ppins.metazoan import MetazoanImporter
from importers.ppins.string import STRINGImporter
import importers.ppins.proteinids as pids
import csv

ppins_dict = dict()

# See also https://github.com/suns-chung/ShortLoopCommonality/tree/master


def new_entry() -> dict:
    entry = {
        "interaction_pair": list(),
        "experiment_type": list(),
        "string_experimental_score": -1.0,
        "string_combined_score": -1.0,
        "metazoan_score": -1.0,
        "bioplex_p_int": -1.0,
        "bioplex_p_wrong": -1.0,
        "biogrid_experimental_system": list(),
        "biogrid_experimental_system_type": list(),
        "biogrid_throughput": list(),
        "biogrid_score": list(),
        "dip_detection_method": list(),
        "publication": list(),
        "sources": set(),
    }
    return entry


class PPINs:
    """Class that encapsulates routines to integrate various protein-protein interaction network datasets"""

    def __init__(
        self,
        node_protein_filename: str,
        ppins_dataset_folder: str,
        rel_interacts_filename: str,
        uniprot_ids_filename: str,
    ):
        self.node_protein_filename = node_protein_filename
        self.rel_interacts_filename = rel_interacts_filename
        self.uniprot_ids_filename = uniprot_ids_filename
        self.biogrid_filename = f"{ppins_dataset_folder}/biogrid/BIOGRID-ORGANISM-Homo_sapiens-4.4.236.tab3.txt"
        self.bioplex_filename = (
            f"{ppins_dataset_folder}/bioplex/BioPlex_interactionList_v4a.tsv"
        )
        self.cell_y2h_hi_union_filename = (
            f"{ppins_dataset_folder}/cell_y2h_hi_union/HI-union.tsv"
        )
        self.dip_filename = f"{ppins_dataset_folder}/dip/Hsapi20170205.txt"
        self.intact_folder = f"{ppins_dataset_folder}/intact/species"
        self.metazoan_filename = f"{ppins_dataset_folder}/metazoan/Hsapiens_projected_PPI_metazoan_conserved_map.txt"
        self.string_filename = (
            f"{ppins_dataset_folder}/string/9606.protein.links.detailed.v12.0.txt"
        )

        # data structures
        self.ppins_dict = dict()
        self.uniprot_primary_secondary_accession_dict = dict()
        self.uniprot_secondary_primary_accession_dict = dict()

    def _load_uniprot_accessions(self):
        """Loads protein accession IDs for filtering purposes"""
        node_protein_csv = ICSV(
            filename=self.node_protein_filename,
            file_type="node_protein",
            max_items=-1,
            delimiter=";",
        )
        # list of chars to be replaced in multi-valued accession_ids
        replace_chars = {"'": "", "[": "", "]": "", " ": ""}
        self.uniprot_primary_secondary_accession_dict = dict()
        self.uniprot_secondary_primary_accession_dict = dict()
        while True:
            item = node_protein_csv.get_next_item()
            if item is None:
                node_protein_csv.close()
                return

            primary_id = item["accession:ID(protein-id)"]
            # TODO: May need refactor due to change in format of list of strings from node_protein.csv
            accession_ids = item["alternative_accessions:string[]"]
            # remove '[] and space
            accession_ids = "".join(replace_chars.get(s, s) for s in accession_ids)
            accession_ids = accession_ids.split(",")
            self.uniprot_primary_secondary_accession_dict[primary_id] = accession_ids
            for secondary_id in accession_ids:
                self.uniprot_secondary_primary_accession_dict[secondary_id] = primary_id

    def create_unified_ppin(self):
        # retrieve the list of uniprot accession IDs (primary and ancillary joined)
        self._load_uniprot_accessions()

        pids.load_uniprot_ids(uniprot_ids_filename=self.uniprot_ids_filename)

        # import biogrid interactions
        biogridimp = BioGRIDImporter(biogrid_filename=self.biogrid_filename)
        biogridimp.merge_entries()

        # import IntAct interactions
        # TODO: Clean up the code of IntActImporter
        intactimp = IntActImporter(intact_folder=self.intact_folder)
        intactimp.merge_entries()

        # import STRING interactions
        stringimp = STRINGImporter(string_filename=self.string_filename)
        stringimp.merge_entries()

        # import DIP
        dipimp = DIPImporter(dip_filename=self.dip_filename)
        dipimp.merge_entries()

        # Large-scale studies (MS experiments)

        # NOTE: hspn (A census of human soluble protein complexes): Havugimana et al. (2012)
        # https://thebiogrid.org/142063/publication/a-census-of-human-soluble-protein-complexes.html
        # is not imported as it has been merged into BioGrid. See {ppins_dataset_folder}/hspn/

        # import Metaozan interactions
        metazoanimp = MetazoanImporter(metazoan_filename=self.metazoan_filename)
        metazoanimp.merge_entries()

        # import cell_y2h_hi_union
        cell_y2h_hi_unionimp = CellY2HHIUnionImporter(
            cell_y2h_hi_union_filename=self.cell_y2h_hi_union_filename
        )
        cell_y2h_hi_unionimp.merge_entries()

        # import BioPlex interactions
        biopleximp = BioPlexImporter(bioplex_filename=self.bioplex_filename)
        biopleximp.merge_entries()

        # NOTE: HeinMS (A Human Interactome in Three Quantitative Dimensions Organized by Stoichiometries and Abundances): Hein et al. (2015)
        # https://thebiogrid.org/188719/publication/a-human-interactome-in-three-quantitative-dimensions-organized-by-stoichiometries-and-abundances.html
        # is not imported as it has been merged into BioGrid. See {ppins_dataset_folder}/heinms/

        # dump merged dataset into the rel_interacts csv file
        with open(self.rel_interacts_filename, "w") as ppi_file:
            ppi_csv = csv.writer(ppi_file, delimiter=";", quotechar="'")
            ppi_csv.writerow(
                [
                    ":START_ID(protein-id)",
                    ":END_ID(protein-id)",
                    "experiment_type:string[]",
                    "string_experimental_score:double",
                    "string_combined_score:double",
                    "metazoan_score:double",
                    "bioplex_p_int:double",
                    "bioplex_p_wrong:double",
                    "biogrid_experimental_system:string[]",
                    "biogrid_experimental_system_type:string[]",
                    "biogrid_throughput:string[]",
                    "biogrid_score:double[]",
                    "dip_detection_method:string[]",
                    "publication:string[]",
                    "sources:string[]",
                    ":TYPE",
                ]
            )
            for entry in ppins_dict.values():
                # this format produces double double quotes
                # experiment_type_str = ",".join(
                #     f"'{w}'" for w in entry["experiment_type"]
                # )
                experiment_type_str = ",".join(entry["experiment_type"])
                biogrid_scores_str = ",".join([str(x) for x in entry["biogrid_score"]])
                biogrid_experimental_system_str = ",".join(
                    entry["biogrid_experimental_system"]
                )
                biogrid_experimental_system_type_str = ",".join(
                    entry["biogrid_experimental_system_type"]
                )
                biogrid_throughput_str = ",".join(entry["biogrid_throughput"])
                dip_detection_method_str = ",".join(entry["dip_detection_method"])
                publication_str = ",".join(entry["publication"])
                sources_str = ",".join(sorted(entry["sources"]))
                ppi_csv.writerow(
                    [
                        entry["interaction_pair"][0],
                        entry["interaction_pair"][1],
                        experiment_type_str,
                        entry["string_experimental_score"],
                        entry["string_combined_score"],
                        entry["metazoan_score"],
                        entry["bioplex_p_int"],
                        entry["bioplex_p_wrong"],
                        biogrid_experimental_system_str,
                        biogrid_experimental_system_type_str,
                        biogrid_throughput_str,
                        biogrid_scores_str,
                        dip_detection_method_str,
                        publication_str,
                        sources_str,
                        "interacts",
                    ]
                )
