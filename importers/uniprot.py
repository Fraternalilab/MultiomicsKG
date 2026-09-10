import time
import datetime
import csv
import xml.etree.ElementTree as ET
from collections import namedtuple

ns = {"ns": "http://uniprot.org/uniprot"}

NS = "{http://uniprot.org/uniprot}"

Tags = namedtuple("Tags", "tag tagns")

entry = Tags("entry", NS + "entry")
name = Tags("name", NS + "name")
accession = Tags("accession", NS + "accession")
organism = Tags("organism", NS + "organism")

skip_items = 0
max_items = -1
estimated_number_items = 20500  # approx. number of items in uniprot_sprot.xml


class UniprotParser:
    """Represents a Uniprot file parser"""

    def __init__(self, filename):
        self.items_read = 0
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

    def get_attribute_value(self, element: ET.Element, path: str, attr: str):
        value = ""
        elements = element.findall(path=path, namespaces=ns)
        if len(elements) > 1:
            raise Exception("More than one AlphaFold ID per protein entry")
        if len(elements) == 1:
            value = elements[0].attrib[attr]
        return value

    def get_organism(self, element: ET.Element):
        item = element.find("ns:organism", namespaces=ns)
        if item is None:
            return None
        else:
            # get only NCBI Taxonomy items
            dbreference = item.find(
                "ns:dbReference[@type='NCBI Taxonomy']", namespaces=ns
            )
            if dbreference is None:
                raise Exception("Protein has no organism in NCBI Taxonomy")
            organism = {}
            organism["scientific_name"] = item.findtext(
                "ns:name[@type='scientific']", default="-", namespaces=ns
            )
            organism["name"] = item.findtext(
                "ns:name[@type='common']", default="-", namespaces=ns
            )
            organism["id"] = dbreference.get("id")
            organism["source"] = "NCBI Taxonomy"
            return organism

    def get_next_item(self) -> namedtuple:
        if self.items_read == max_items:
            return None

        for event, element in self.ittree:
            # print("{event} : {tag} -> {element}".format(event=event, tag=element.tag, element=element))
            if element.tag == entry.tagns:
                self.items_read += 1

                # if self.items_read % 100 == 0:
                #     print(
                #         "Uniprot items read: {items_read}".format(
                #             items_read=self.items_read
                #         ),
                #         end="\r",
                #     )

                if self.items_read < skip_items:
                    element.clear()
                    continue

                UniprotItem = namedtuple("UniprotItem", "protein organism model")
                protein = {}

                # ns:accession contains all accession IDs, take the 1st one as principal
                ns_accession = self.get_text(element=element, path="ns:accession")
                accession_ids = ns_accession.split(",")

                protein["accession"] = accession_ids[0]
                protein["alternative_accessions"] = accession_ids

                protein["name"] = self.get_text(element=element, path="ns:name")
                protein["full_name"] = self.get_text(
                    element=element, path="ns:protein/ns:recommendedName/ns:fullName"
                )

                organism = self.get_organism(element=element)

                sequence_el = element.find("ns:sequence", namespaces=ns)
                protein["sequence_length"] = sequence_el.get("length")
                protein["sequence_mass"] = sequence_el.get("mass")
                protein["sequence"] = sequence_el.text

                protein["function"] = self.get_text(
                    element=element, path='ns:comment[@type="function"]/ns:text'
                )

                model_id = self.get_attribute_value(
                    element=element,
                    path='ns:dbReference[@type="AlphaFoldDB"]',
                    attr="id",
                )

                model = {"id": model_id}

                # if type(protein["name"]) == list:
                #     print("MULTIPLE NAMES!!!!")
                # if type(protein["accession"]) == list:
                #     print("MULTIPLE ACCESSION!!!")

                # print("Acc: {acc} | Name: {name}".format(acc=elAccession[0].text, name=elName[0].text))
                element.clear()
                return UniprotItem(protein, organism, model)

        return None

    def close(self):
        print(
            "\n{module} read {numitems} items".format(
                module=__package__, numitems=self.items_read
            )
        )


class UniprotImporter:
    """
    Represents a Uniprot importer
    https://www.uniprot.org/help/downloads
    """

    def __init__(
        self,
        filename,
        node_protein_filename,
        node_organism_filename,
        rel_sampled_from_filename,
        node_model_filename,
        rel_protein_model_filename,
        human_only,
        estimated_number_items=estimated_number_items
    ):
        self.data = []
        self.parser = UniprotParser(filename)
        self.node_protein_filename = node_protein_filename
        self.node_organism_filename = node_organism_filename
        self.rel_sampled_from_filename = rel_sampled_from_filename
        self.node_model_filename = node_model_filename
        self.rel_protein_model_filename = rel_protein_model_filename
        self.estimated_number_items = estimated_number_items
        self.human_only = human_only
        self.num_imported = 0
        self.report_line = (
            "Uniprot items imported: {num_items:7d}   |   "
            + "Items/s: {speed:8.2f}   |   "
            + "Time elapsed: {time_elapsed}   |   "
            + "ETA: {eta}"
        )

    def close(self):
        self.parser.close()

    def importdb(self):
        start_counter = time.perf_counter()
        stop_counter = start_counter
        start = time.time()

        organism_keys = set()

        model_keys = set()

        node_protein_file = open(self.node_protein_filename, "w")
        node_organism_file = open(self.node_organism_filename, "w")
        rel_sampled_from_file = open(self.rel_sampled_from_filename, "w")
        node_model_file = open(self.node_model_filename, "w")
        rel_protein_model_file = open(self.rel_protein_model_filename, "w")
        node_protein_csv = csv.writer(node_protein_file, delimiter=";", quotechar="'")
        node_organism_csv = csv.writer(node_organism_file, delimiter=";", quotechar="'")
        rel_sampled_from_csv = csv.writer(
            rel_sampled_from_file, delimiter=";", quotechar="'"
        )
        node_model_csv = csv.writer(node_model_file, delimiter=";", quotechar="'")
        rel_protein_model_csv = csv.writer(
            rel_protein_model_file, delimiter=";", quotechar="'"
        )

        node_protein_csv.writerow(
            [
                "accession:ID(protein-id)",
                "alternative_accessions:string[]",
                "name:string",
                "full_name:string",
                "sequence_length:long",
                "sequence_mass:long",
                "sequence:string",
                "function:string",
            ]
        )
        node_organism_csv.writerow(
            ["scientific_name", "name", "id:ID(organism-id)", "source"]
        )
        rel_sampled_from_csv.writerow(
            [":START_ID(protein-id)", ":END_ID(organism-id)", ":TYPE"]
        )

        node_model_csv.writerow(["id:ID(model-id)"])
        rel_protein_model_csv.writerow(
            [
                ":START_ID(protein-id)",
                "chain",
                "res_beg",
                "res_end",
                "pdb_beg",
                "pdb_end",
                "sp_beg",
                "sp_end",
                ":END_ID(model-id)",
                ":TYPE",
            ]
        )

        while True:
            item = self.parser.get_next_item()

            if item is None:
                print("\nNo more items")
                node_protein_file.close()
                node_organism_file.close()
                rel_sampled_from_file.close()
                node_model_file.close()
                rel_protein_model_file.close()
                break

            # check if source organism is human
            if self.human_only and item.organism["id"] != "9606":
                continue

            # check for existing protein
            # accession = item.protein["accession"]
            # proteinvx = self.persistence.get_vertex(label="protein", key_name="accession", key_value=accession)
            # if len(proteinvx) == 0:
            #     # create a new protein
            # proteinvx = self.persistence.create_vertex(label="protein", props=item.protein)

            # check for existing organism
            # organism_id = item.organism["id"]
            # organismvx = self.persistence.get_vertex(label="organism", key_name="id", key_value=organism_id)
            # if len(organismvx) == 0:
            # create organism
            # organismvx = self.persistence.create_vertex(label="organism", props=item.organism)

            # create protein, organism and link them if pattern does not exist yet

            # self.persistence.create_linked_pair(
            #     labelA="protein",
            #     keyA="accession",
            #     itemA=item.protein,
            #     labelB="organism",
            #     keyB="id",
            #     itemB=item.organism,
            #     relation="sampled_from",
            # )
            node_protein_csv.writerow(
                [
                    item.protein["accession"],
                    ",".join(item.protein["alternative_accessions"]),
                    item.protein["name"],
                    item.protein["full_name"],
                    item.protein["sequence_length"],
                    item.protein["sequence_mass"],
                    item.protein["sequence"],
                    item.protein["function"],
                ]
            )
            if item.organism["id"] not in organism_keys:
                node_organism_csv.writerow(item.organism.values())
                organism_keys.add(item.organism["id"])
            rel_sampled_from_csv.writerow(
                [item.protein["accession"], item.organism["id"], "sampled_from"]
            )

            model_id = item.model["id"]
            if model_id != "":
                uniprot_id = item.protein["accession"]

                map_protein_pdb = dict()

                map_protein_pdb["uniprot"] = uniprot_id
                map_protein_pdb["chain"] = "A"
                map_protein_pdb["res_beg"] = 1
                map_protein_pdb["res_end"] = int(item.protein["sequence_length"])
                map_protein_pdb["pdb_beg"] = 1
                map_protein_pdb["pdb_end"] = int(item.protein["sequence_length"])
                map_protein_pdb["sp_beg"] = 1
                map_protein_pdb["sp_end"] = int(item.protein["sequence_length"])
                map_protein_pdb["model_id"] = model_id
                map_protein_pdb["type"] = "maps"

                # PDB node
                if model_id not in model_keys:
                    node_model_csv.writerow([model_id])
                    model_keys.add(model_id)

                # PROTEIN->PDB relation
                rel_protein_model_csv.writerow(map_protein_pdb.values())

            self.num_imported += 1
            if self.num_imported % 1000 == 0:
                stop_counter = time.perf_counter()
                speed = 1000 / (stop_counter - start_counter)
                eta = str(
                    datetime.timedelta(
                        seconds=round(self.estimated_number_items / speed)
                    )
                )
                time_elapsed = str(
                    datetime.timedelta(seconds=round(time.time() - start))
                )
                print(
                    self.report_line.format(
                        num_items=self.num_imported,
                        speed=speed,
                        time_elapsed=time_elapsed,
                        eta=eta,
                    ),
                    end="\r",
                )
                start_counter = time.perf_counter()
