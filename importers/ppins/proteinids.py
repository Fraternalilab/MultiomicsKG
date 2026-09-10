import csv
import sys

# contains set of cross-reference IDs of proteins by Uniprot Accession
uniprot_ids_dict = dict()
# contains name of proteins (Uniprot ID) by Uniprot Accession
uniprot_names_dict = dict()
# contains Uniprot Accession by names of proteins
uniprot_acc_by_names_dict = dict()
# contains Ensembl Gene ids by Uniprot Accession
ensembl_gen_dict = dict()
# contains Ensembl Protein ids by Uniprot Accession
ensembl_pro_dict = dict()

# Uniprot ID Mapping
# https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/by_organism/HUMAN_9606_idmapping_selected.tab.gz


class ProteinIDSet:

    def __init__(
        self,
        uniprot_ac: str,
        uniprot_id: str,
        entrez_gene: str,
        pdb: str,
        ensembl_gen: str,
        ensembl_pro: str,
        ensembl_tra: str,
        go_ids: str,
    ):
        self.uniprot_ac = uniprot_ac
        self.uniprot_id = uniprot_id
        self.entrez_gene = entrez_gene
        self.ensembl_gen = ensembl_gen
        self.ensembl_pro = ensembl_pro
        self.ensembl_tra = ensembl_tra
        self.go_ids = go_ids


def load_uniprot_ids(uniprot_ids_filename: str):

    csv.field_size_limit(sys.maxsize)
    with open(uniprot_ids_filename) as csvfile:
        uniprot_ids_csv = csv.reader(csvfile, delimiter="\t")
        for row in uniprot_ids_csv:
            # see header description here: https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/README
            protein_id_set = ProteinIDSet(
                uniprot_ac=row[0],
                uniprot_id=row[1],
                entrez_gene=row[2],
                pdb=row[5],
                ensembl_gen=row[18],
                ensembl_pro=row[20],
                ensembl_tra=row[19],
                go_ids=row[6],
            )
            uniprot_names_dict[protein_id_set.uniprot_ac] = protein_id_set.uniprot_id
            uniprot_acc_by_names_dict[protein_id_set.uniprot_id] = (
                protein_id_set.uniprot_ac
            )
            uniprot_ids_dict[protein_id_set.uniprot_ac] = protein_id_set
            # buildEntrezUniprotIDLib(
            #     aProt.getEntrezID(), aProt.getUniprotAC(), aProt.getUniprotID()
            # )
            # buildENSTUniprotLib(
            #     aProt.getEnsemblTraID(), aProt.getUniprotAC(), aProt.getUniprotID()
            # )
            create_ensembl_gen_to_uniprotacc(
                protein_id_set.ensembl_gen,
                protein_id_set.uniprot_ac,
                protein_id_set.uniprot_id,
            )
            create_ensembl_pro_to_uniprotacc(
                protein_id_set.ensembl_pro,
                protein_id_set.uniprot_ac,
                protein_id_set.uniprot_id,
            )


def create_ensembl_gen_to_uniprotacc(ensg_ids, uniprot_acc, uniprot_id):
    for id in ensg_ids.split(";"):
        ensg_id_with_suffix = id.strip()
        ensg_id = ensg_id_with_suffix.split(".")[0]
        if len(ensg_id) > 0:
            if ensg_id not in ensembl_gen_dict.keys():
                ensembl_gen_dict[ensg_id] = uniprot_acc
            else:
                if not uniprot_id.startswith(uniprot_acc):
                    ensembl_gen_dict[ensg_id] = uniprot_acc


def create_ensembl_pro_to_uniprotacc(ensp_ids, uniprot_acc, uniprot_id):
    for id in ensp_ids.split(";"):
        ensp_id_with_suffix = id.strip()
        ensp_id = ensp_id_with_suffix.split(".")[0]
        if len(ensp_id) > 0:
            if ensp_id not in ensembl_pro_dict.keys():
                ensembl_pro_dict[ensp_id] = uniprot_acc
            else:
                if not uniprot_id.startswith(uniprot_acc):
                    ensembl_pro_dict[ensp_id] = uniprot_acc
        # else:
        #     print(
        #         f"Uniprot accession {uniprot_acc} does not contain a mapped Ensembl Protein ID"
        #     )


def findUniprotID(uacc):
    if uacc in uniprot_ids_dict.keys():
        return uniprot_ids_dict[uacc].uniprot_id
    else:
        return "NA"


def confirmUniprotID(uid):
    if not uid.endswith("_HUMAN"):
        uid_cols = uid.split("-")
        # 		print(uid_cols[0])
        return findUniprotID(uid_cols[0].upper())
    else:
        if uid in uniprot_names_dict.values():
            return uid
        else:
            return "NA"
