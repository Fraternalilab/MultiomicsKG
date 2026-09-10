import logging
import logging.handlers
import argparse
import filenames
from procedures.structural_tools import StructuralToolsRunner

if __name__ == "__main__":
    logging.basicConfig(
        filename="sdgenerator.log",
        filemode="w",
        format="%(asctime)s %(processName)-10s %(name)s %(levelname)-8s %(message)s",
        level=logging.DEBUG,
    )

    parser = argparse.ArgumentParser(
        description="Process PDB files and generate data from structural analysis tools",
        argument_default=argparse.SUPPRESS,
    )
    parser.add_argument(
        "-d",
        required=True,
        type=str,
        dest="datastore_prefix",
        help="path of the data store, e.g.: /Users/myuser/datastore",
    )
    parser.add_argument(
        "-s",
        required=True,
        type=str,
        dest="scrap_root",
        help="path of the root for the scrap folder, e.g.: /Users/myuser/other_drive",
    )
    parser.add_argument(
        "-p",
        required=False,
        type=str,
        default="",
        dest="run_prefix",
        help='prefix to label this run, e.g.: "test_"',
    )
    parser.add_argument(
        "-m",
        required=True,
        type=int,
        default=0,
        dest="homodimers",
        help="run the process targeting the homodimers dataset, e.g.: 0 (false)",
    )
    args = parser.parse_args()

    fn = filenames.FileNames(args.datastore_prefix, args.scrap_root, args.run_prefix)
    homodimers = bool(args.homodimers)

    pdb_structural_tools_runner = StructuralToolsRunner(
        node_pdb_filename=fn.node_pdb_filename,
        pdb_folder=fn.pdb_database_root,
        structure_db_type="biounit",
        scrap_folder=fn.scrap_folder,
        issues_filename=fn.pdb_issues_filename,
        pdb_residues_progress_filename=fn.pdb_residues_progress_filename,
        skip_defects=True,
        homodimers=False,
        use_cache=True,
    )
    pdb_structural_tools_runner.run()

    model_structural_tools_runner = StructuralToolsRunner(
        node_pdb_filename=fn.node_model_filename,
        pdb_folder=fn.alphafold_database_root,
        structure_db_type="alphafold",
        scrap_folder=fn.scrap_folder,
        issues_filename=fn.model_issues_filename,
        pdb_residues_progress_filename=fn.model_residues_progress_filename,
        skip_defects=True,
        homodimers=False,
        use_cache=True,
    )
    model_structural_tools_runner.run()

    # # run node_pdb to node_residue layer importer
    # # reads: node_hdpdb.csv, homodimers pdb dataset
    # # writes: node_chain.csv, node_residue.csv, rel_pdb_chain.csv, rel_pdb_residue.csv,
    # #         rel_chain_residue.csv, issues.json, physio.json, nonphysio.json
    # # pdb_residues = PDBResidues(
    # #     node_hdpdb_datafile, homodimers_pdb_database_root, scrap_folder
    # # )
    # pdb_residues = StructureProcessor(
    #     node_pdb_filename=fn.node_pdb_filename,
    #     pdb_folder=fn.pdb_database_root,
    #     scrap_folder=fn.scrap_folder,
    #     issues_filename=fn.issues_filename,
    #     physio_filename=fn.physio_filename,
    #     nonphysio_filename=fn.nonphysio_filename,
    #     node_chain_filename=fn.node_chain_filename,
    #     node_residue_filename=fn.node_residue_filename,
    #     rel_pdb_chain_filename=fn.rel_pdb_chain_filename,
    #     rel_pdb_residue_filename=fn.rel_pdb_residue_filename,
    #     rel_chain_residue_filename=fn.rel_chain_residue_filename,
    #     pdb_residues_progress_filename=fn.pdb_residues_progress_filename,
    #     homodimers=homodimers,
    # )
    # pdb_residues.run()
    # pdb_residues.close()
