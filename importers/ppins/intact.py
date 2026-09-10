import xml.etree.ElementTree as ET
import glob
import importers.ppins.proteinids as pids
import importers.ppins.ppins as ppins


class IntActImporter:
    """
    Represents an IntAct importer
    https://www.ebi.ac.uk/intact/home
    """

    def __init__(self, intact_folder: str):
        self.intact_folder = intact_folder
        self.protIDmap = dict()

    def importfile(self):
        pass

    def merge_entries(self):

        # TODO: which files to consider?
        xmlfiles = glob.glob(f"{self.intact_folder}/[Hh]uman*[0-9].xml")
        for afile in xmlfiles:
            print(afile)
            self._buildProtPsimi(afile)
            self._buildPPIPsimi(afile)

    def _buildProtPsimi(self, file_name: str):
        ns = {"psimins": "http://psi.hupo.org/mi/mif"}
        tree = ET.parse(file_name)
        root = tree.getroot()
        for interactor in root.findall(".//psimins:interactor", ns):
            intid = interactor.get("id")
            names = interactor.find("psimins:names", ns)
            # 	print(intid)
            label = names.find("psimins:shortLabel", ns).text
            # UniprotID validation
            self.protIDmap[intid] = pids.confirmUniprotID(label.upper())
        # 	print(label)

    def _buildPPIPsimi(self, file_name: str):
        # 	import xml.etree.ElementTree as ET
        ns = {"psimins": "http://psi.hupo.org/mi/mif"}
        tree = ET.parse(file_name)
        root = tree.getroot()
        for interactor in root.findall(".//psimins:interactor", ns):
            names = interactor.find("psimins:names", ns)
            # 		print(names.tag)
            label = names.find("psimins:shortLabel", ns).text
        # 		print(label)
        for interaction in root.findall(".//psimins:interaction", ns):
            names = interaction.find("psimins:names", ns)
            label = names.find("psimins:shortLabel", ns).text
            # 		print(label)
            # 		for child in interaction :
            # 			print(child.tag)
            ppi_pair = []
            for participants in interaction.findall(".//psimins:participant", ns):
                protname = participants.find("psimins:names", ns)
                if protname is not None:
                    protlabel = protname.find("psimins:shortLabel", ns).text
                    # 				print(label)
                    interactorRef = participants.find("psimins:interactorRef", ns).text
                    tuid = self.protIDmap[interactorRef]
                    # 				if not tuid.endswith('_HUMAN') :
                    # 					tuid_cols=tuid.split('-')
                    # 					uid=findUniprotID(tuid_cols[0].upper())
                    # 	if uid.endswith('_HUMAN') :
                    # 					print(tuid+'   '+uid)
                    ppi_pair.append(tuid)
                    # 	else :
                    # 		ppi_pair.append(tuid)
            # 				else :
            # 					uid=confirmUniprotID(tuid.upper())
            # 					print(tuid+'   '+uid)
            # 					ppi_pair.append(uid)
            if len(ppi_pair) >= 2:
                if ppi_pair[0].endswith("_HUMAN"):
                    for i in range(1, len(ppi_pair)):
                        # 					if frozenset(ppi_pair[0],ppi_pair[i]) not in IntActPPIpool :
                        # Ignore Self-loop
                        if ppi_pair[i].endswith("_HUMAN"):
                            protein1_uniacc = pids.uniprot_acc_by_names_dict.get(
                                ppi_pair[0], None
                            )
                            protein2_uniacc = pids.uniprot_acc_by_names_dict.get(
                                ppi_pair[i], None
                            )
                            if (protein1_uniacc is not None) and (
                                protein2_uniacc is not None
                            ):
                                if (
                                    protein1_uniacc not in pids.uniprot_ids_dict
                                    or protein2_uniacc not in pids.uniprot_ids_dict
                                ):
                                    print(
                                        f"Unexpected interaction bewtween {protein1_uniacc} and {protein2_uniacc}"
                                    )
                                else:
                                    interaction_pair = sorted(
                                        [protein1_uniacc, protein2_uniacc]
                                    )
                                    key = f"{interaction_pair[0]}_{interaction_pair[1]}"

                                    entry = dict()
                                    if key in ppins.ppins_dict:
                                        entry = ppins.ppins_dict[key]
                                    else:
                                        entry = ppins.new_entry()
                                    entry["interaction_pair"] = interaction_pair
                                    entry["sources"].add("intact")
                                    ppins.ppins_dict[key] = entry
                            # else:
                            #     print(f"Item missing at least one interactor's Uniprot Accession: {item}")
