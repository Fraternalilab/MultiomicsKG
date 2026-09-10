from parsers.icsv import ICSV


class HeinMSImporter(ICSV):
    """
    Represents a HeinMS dataset importer
    https://thebiogrid.org/188719/publication/a-human-interactome-in-three-quantitative-dimensions-organized-by-stoichiometries-and-abundances.html

    NOTE: Not used as already imported into BioGrid
    """

    def __init__(self, string_filename: str):
        super().__init__(
            filename=string_filename,
            file_type="heinms",
            max_items=-1,
            delimiter="\t",
            header_filter_prefix="",
            compressed=False,
        )
        self.entries = list()
        self.entrez_pool = list()

    def importfile(self):
        pass

    def get_entries(self) -> list:
        self.entries = list()
        # performs standard import (calls process_item())
        super().importfile()
        super().close()
        return self.entries

    def process_item(self, item: dict):

        if item[0].startswith("entrez gene/locuslink:") and item[1].startswith(
            "entrez gene/locuslink:"
        ):
            protein1 = item[0].split(":")
            protein2 = item[1].split(":")
        else:
            print(item[0] + " " + item[1])

        if len(protein1) >= 2 and len(protein2) >= 2:
            # 	print(prot1[1] + ' ' + prot2[1])

            protein1_uniprot_id = getUniprotIDEntrez(protein1[1].strip())
            protein2_uniprot_id = getUniprotIDEntrez(protein2[1].strip())

            if (protein1_uniprot_id != "NA") and (protein2_uniprot_id != "NA"):
                entry = {
                    "protein1": protein1_uniprot_id,
                    "protein2": protein2_uniprot_id,
                }
                self.entries.append(entry)
            else:
                entry = {
                    "protein1": protein1_uniprot_id,
                    "protein2": protein2_uniprot_id,
                }
                self.entrez_pool.append(entry)
