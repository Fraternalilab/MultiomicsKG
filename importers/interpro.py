from parsers.icsv import ICSV
import xml.etree.ElementTree as ET
from collections import namedtuple
import csv

ns = {"ns": ""}

NS = ""

Tags = namedtuple("Tags", "tag tagns")

interpro = Tags("interpro", NS + "interpro")
name = Tags("name", NS + "name")


class InterproParser:
    """Represents a Uniprot file parser"""

    def __init__(
        self, filename: str, max_items: int = -1, skip_items: int = -1
    ) -> None:
        self.items_read = 0
        self.max_items = max_items
        self.skip_items = skip_items
        self.ittree = ET.iterparse(filename)

    def get_text(self, element: ET.Element, path: str):
        values = element.findall(path, namespaces=ns)
        text = ""
        i = 0
        j = len(values)
        while i < j:
            text = text + values[i].text
            i += 1
            if i < j:
                text = text + ","

        return text

    def importdb(self) -> dict:
        interpro_items = dict()
        if self.items_read == self.max_items:
            return interpro_items

        for _, element in self.ittree:
            if element.tag == interpro.tagns:
                self.items_read += 1

                if self.items_read < self.skip_items:
                    element.clear()
                    continue

                interpro_type = element.attrib["type"]
                # import only Domain, Family and Homologous_superfamily types
                if interpro_type in [
                    "Domain",
                    "Family",
                    "Homologous_superfamily",
                ]:
                    interpro_id = element.attrib["id"]
                    child_list = element.findall("child_list/rel_ref")
                    children = []
                    for child in child_list:
                        children.append(child.attrib["ipr_ref"])
                    interpro_item = {
                        "short_name": element.attrib["short_name"],
                        "name": self.get_text(element=element, path="name"),
                        "type": interpro_type,
                        "is-llm": element.attrib["is-llm"] == "true",
                        "is-llm-reviewed": element.attrib["is-llm-reviewed"] == "true",
                        "children": children,
                    }

                    interpro_items[interpro_id] = interpro_item
            if self.items_read % 1000 == 0:
                print(
                    "interpro items read: {num_items:7d}".format(
                        num_items=self.items_read
                    ),
                    end="\r",
                )
        return interpro_items

    def close(self):
        print(
            "\n{module} read {numitems} items".format(
                module=__package__, numitems=self.items_read
            )
        )


class InterproImporter:
    """
    Represents an Interpro importer
    https://www.ebi.ac.uk/interpro/download/
    """

    def __init__(
        self,
        protein2interpro_filename: str,
        interpro_filename: str,
        entrylist_filename: str,
        node_interpro_entries_filename: str,
        rel_interpro_filename: str,
        rel_protein_interpro_entries_filename: str,
        node_protein_filename: str,
    ) -> None:
        self.data = []
        self.protein2interpro_filename = protein2interpro_filename
        self.interpro_filename = interpro_filename
        self.entrylist_filename = entrylist_filename
        self.node_interpro_entries_filename = node_interpro_entries_filename
        self.rel_interpro_filename = rel_interpro_filename
        self.rel_protein_interpro_entries_filename = (
            rel_protein_interpro_entries_filename
        )
        self.node_protein_filename = node_protein_filename
        self.uniprot_primary_secondary_accession_dict = dict()
        self.uniprot_secondary_primary_accession_dict = dict()
        self.protein_interpro_map = dict()
        self.interpro_entry_map = dict()
        self.interpro_map = dict()

    def _load_uniprot_accessions(self):
        """Loads protein accession IDs for filtering purposes"""
        node_protein_csv = ICSV(
            filename=self.node_protein_filename,
            file_type="node_protein",
            max_items=-1,
            delimiter=";"
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

    def _load_protein_to_interpro(self) -> None:
        # Load Protein2Interpro entries
        protein2interpro_csv = ICSV(
            filename=self.protein2interpro_filename,
            file_type="protein2interpro",
            delimiter="\t",
            compressed=True,
            skip_rows=0,
            fieldnames=[
                "protein_accession",
                "interpro_id",
                "name",
                "source_id",
                "start",
                "end",
            ],
        )

        num_interpro_entries = 0
        while True:
            item = protein2interpro_csv.get_next_item()
            if item is None:
                protein2interpro_csv.close()
                break
            else:
                # check if protein accession is in uniprot accessions
                if (
                    item["protein_accession"]
                    in self.uniprot_secondary_primary_accession_dict
                ):
                    num_interpro_entries += 1
                    primary_accession = self.uniprot_secondary_primary_accession_dict[
                        item["protein_accession"]
                    ]
                    interpro_item = {
                        "interpro_id": item["interpro_id"],
                        "name": item["name"],
                        "source_id": item["source_id"],
                        "start": int(item["start"]),
                        "end": int(item["end"]),
                    }
                    if primary_accession in self.protein_interpro_map:
                        self.protein_interpro_map[primary_accession].append(
                            interpro_item
                        )
                    else:
                        self.protein_interpro_map[primary_accession] = []
                        self.protein_interpro_map[primary_accession].append(
                            interpro_item
                        )

        print(
            f"Imported {num_interpro_entries} InterPro entries for {len(self.protein_interpro_map)} proteins."
        )

    def _load_interpro(self) -> None:
        """Loads Interpro entries from interpro.xml.gz"""
        interpro_parser = InterproParser(filename=self.interpro_filename)
        self.interpro_map = interpro_parser.importdb()
        interpro_parser.close()
        print(f"Imported {len(self.interpro_map)} InterPro entries from XML file.")

    def _load_interpro_entrylist(self) -> None:
        """Loads Interpro entry list"""
        interpro_entry_csv = ICSV(
            filename=self.entrylist_filename,
            file_type="interpro_entrylist",
            max_items=-1,
            delimiter="\t",
            skip_rows=1,
            fieldnames=[
                "ENTRY_AC",
                "ENTRY_TYPE",
                "ENTRY_NAME",
            ],
        )
        while True:
            item = interpro_entry_csv.get_next_item()
            if item is None:
                interpro_entry_csv.close()
                break
            # import only Domain, Family and Homologous_superfamily types
            if item["ENTRY_TYPE"] in [
                "Domain",
                "Family",
                "Homologous_superfamily",
            ]:
                self.interpro_entry_map[item["ENTRY_AC"]] = {
                    "entry_type": item["ENTRY_TYPE"],
                    "name": item["ENTRY_NAME"],
                }
        print(
            f"Imported detailed records of {len(self.interpro_entry_map)} InterPro entries"
        )

    def importdb(self) -> None:
        """Imports Interpro data from given files"""

        # retrieve the list of uniprot accession IDs (primary and ancillary joined)
        self._load_uniprot_accessions()

        # retrieve mapping of protein to Interpro entries
        self._load_protein_to_interpro()

        # retrieve Interpro data from interpro.xml.gz
        self._load_interpro()

        # # retrieve Interpro entry details
        # self._load_interpro_entrylist()

        node_interpro_entries_file = open(self.node_interpro_entries_filename, "w")
        node_interpro_entries_csv = csv.writer(
            node_interpro_entries_file, delimiter=";", quotechar="'"
        )
        node_interpro_entries_csv.writerow(
            [
                "id:ID(interpro-id)",
                "short_name:string",
                "name:string",
                "type:string",
                "is-llm",
                "is-llm-reviewed",
            ]
        )

        rel_interpro_file = open(self.rel_interpro_filename, "w")
        rel_interpro_csv = csv.writer(rel_interpro_file, delimiter=";", quotechar="'")
        rel_interpro_csv.writerow(
            [":START_ID(interpro-id)", ":END_ID(interpro-id)", ":TYPE"]
        )

        rel_protein_interpro_entries_file = open(
            self.rel_protein_interpro_entries_filename, "w"
        )
        rel_protein_interpro_entries_csv = csv.writer(
            rel_protein_interpro_entries_file, delimiter=";", quotechar="'"
        )
        rel_protein_interpro_entries_csv.writerow(
            [
                ":START_ID(protein-id)",
                ":END_ID(interpro-id)",
                "source_id:string",
                "start:long",
                "end:long",
                ":TYPE",
            ]
        )

        # write protein to interpro entries relationships and save relevant interpro entries
        selected_interpro_entries = dict()
        for protein_accession in self.protein_interpro_map:
            interpro_items = self.protein_interpro_map[protein_accession]
            for interpro_item in interpro_items:
                interpro_id = interpro_item["interpro_id"]
                if interpro_id in self.interpro_map:
                    interpro_data = self.interpro_map[interpro_id]

                    if interpro_id not in selected_interpro_entries:
                        selected_interpro_entries[interpro_id] = interpro_data

                    # write rel_protein_interpro_entries_csv
                    rel_protein_interpro_entries_csv.writerow(
                        [
                            protein_accession,
                            interpro_id,
                            interpro_item["source_id"],
                            interpro_item["start"],
                            interpro_item["end"],
                            "has_interpro_entry",
                        ]
                    )

                    # # write node_interpro_entries.csv
                    # alternative_accessions = []
                    # if interpro_id in self.interpro_entry_map:
                    #     alternative_accessions.append(
                    #         self.interpro_entry_map[interpro_id]["name"]
                    #     )
                    # node_interpro_entries_csv.writerow(
                    #     [
                    #         interpro_id,
                    #         ",".join(alternative_accessions),
                    #         interpro_data["short_name"],
                    #         interpro_data["name"],
                    #         "",  # sequence_length
                    #         "",  # sequence_mass
                    #         "",  # sequence
                    #         "",  # function
                    #     ]
                    # )

        # write interpro to interpro relationships
        for interpro_id in selected_interpro_entries:
            interpro_data = selected_interpro_entries[interpro_id]
            for child_id in interpro_data["children"]:
                if child_id in selected_interpro_entries:
                    rel_interpro_csv.writerow([interpro_id, child_id, "has_child"])

        # write interpro entries
        for interpro_id in selected_interpro_entries:
            interpro_data = selected_interpro_entries[interpro_id]
            # write node_interpro_entries.csv
            node_interpro_entries_csv.writerow(
                [
                    interpro_id,
                    interpro_data["short_name"],
                    interpro_data["name"],
                    interpro_data["type"],
                    interpro_data["is-llm"],
                    interpro_data["is-llm-reviewed"],
                ]
            )
        node_interpro_entries_file.close()
        rel_interpro_file.close()
        rel_protein_interpro_entries_file.close()
