from Bio.PDB import PDBParser
from Bio.PDB import NeighborSearch
from Bio.PDB import Selection
from Bio import BiopythonWarning
import warnings

# Reference:
# https://biopython.org/docs/latest/api/Bio.PDB.NeighborSearch.html
# https://biopython.org/wiki/The_Biopython_Structural_Bioinformatics_FAQ#the-structure-object
# https://www.biostars.org/p/269579/


def generate_contact_residue_list(
    pdb_id: str, pdb_file_path: str, contact_map_cutoff: float
):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", BiopythonWarning)

        contacts_list = list()

        parser = PDBParser()
        structure = parser.get_structure(pdb_id, pdb_file_path)

        # structure -> model(s) -> chain(s) -> residue(s) -> atom(s)

        # take first model (no NMR support here)
        model = structure[0]

        chain_list = Selection.unfold_entities(model, "C")
        for chain in chain_list:
            chain_id = chain.id
            # get only atoms from chain as we're interested only
            # in intra-chain contacts for now
            atom_list = Selection.unfold_entities(chain, "A")
            # filter by only CA atoms
            ca_atom_list = [ca_atom for ca_atom in atom_list if "CA" == ca_atom.id]

            if len(ca_atom_list) == 0:
                # the chain does not contain Ca atoms, so skip it
                continue

            ns = NeighborSearch(ca_atom_list)
            contact_residues = ns.search_all(radius=contact_map_cutoff, level="R")
            for residue_pair in contact_residues:
                if len(residue_pair) > 2:
                    raise Exception(
                        f"{pdb_id}: unexpected number of contact pair elements"
                    )
                res1 = f"{pdb_id}_{chain_id}_{residue_pair[0].id[1]}"
                res2 = f"{pdb_id}_{chain_id}_{residue_pair[1].id[1]}"
                contacts_list.append((res1, res2))

        # contacts_list is a tuple of res1, res2 (strings with residue ids)
        return contacts_list


if __name__ == "__main__":
    cmap = generate_contact_residue_list(
        "1a1x", "/Volumes/PhD/datastore/cmap_scrap/1a1x.pdbs", contact_map_cutoff=7.0
    )
