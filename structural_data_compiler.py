import logging
import argparse
import filenames
from procedures.structural_compiler import StructuralCompilerRunner
from os import makedirs

if __name__ == "__main__":
    logging.getLogger().setLevel(logging.INFO)

    parser = argparse.ArgumentParser(
        description="Compiles data from structural analysis tools into CSV files ready to import into a DBMS",
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
    args = parser.parse_args()

    fn = filenames.FileNames(args.datastore_prefix, args.scrap_root, args.run_prefix)

    makedirs(fn.graph_output, exist_ok=True)

    # generates node_chain.csv, node_residue.csv, rel_pdb_chain.csv, rel_model_chain.csv,
    # rel_pdb_residue.csv, rel_model_residue.csv, rel_chain_residue.csv, rel_residue_residue.csv
    structural_compiler_runner = StructuralCompilerRunner(
        node_pdb_filename=fn.node_pdb_filename,
        node_model_filename=fn.node_model_filename,
        scrap_folder=fn.scrap_folder,
        pdb_issues_filename=fn.pdb_issues_filename,
        model_issues_filename=fn.model_issues_filename,
        node_chain_filename=fn.node_chain_filename,
        node_residue_filename=fn.node_residue_filename,
        rel_pdb_chain_filename=fn.rel_pdb_chain_filename,
        rel_model_chain_filename=fn.rel_model_chain_filename,
        rel_pdb_residue_filename=fn.rel_pdb_residue_filename,
        rel_model_residue_filename=fn.rel_model_residue_filename,
        rel_chain_residue_filename=fn.rel_chain_residue_filename,
        rel_residue_residue_filename=fn.rel_residue_residue_filename,
        skip_defects=True,
    )
    structural_compiler_runner.run()
