import logging
import rpy2.robjects as robjects
import rpy2.rinterface_lib.callbacks as r_cbs
import subprocess
import os
import shutil
import csv
import xmlrpc.client as xmlrpclib
import argparse

# Redirect console output from R to generic logging
# Ref: https://rpy2.github.io/doc/latest/html/callbacks.html

logger = logging.getLogger()


def r_consolewrite_print_cb(s: str):
    logger.info(s)


def r_consolewrite_warnerror_cb(s: str):
    logger.warning(s)


r_cbs.consolewrite_print = r_consolewrite_print_cb
r_cbs.consolewrite_warnerror = r_consolewrite_warnerror_cb


class PDBEncode:
    """Class that interacts with R to run PDBEncode"""

    def __init__(self, scrap_folder: str):
        self.scrap_folder = scrap_folder
        self.pdbe = robjects.r
        self.pdbe.source("procedures/pdbencode_blocks.R")

    def load_structure(self, file_path: str):
        self.structure = self.pdbe.load_structure(file_path=file_path)

    def encode_chain(self, chain: str) -> str:
        output = self.pdbe.encode_chain(chain=chain)
        return output[0]  # discard the encapsulating R object and return a Python str

    def encode_range(self, chain: str, range: list) -> str:
        range_vector = robjects.vectors.IntVector(range)
        output = self.pdbe.encode_range(chain=chain, range=range_vector)
        return output[0]  # discard the encapsulating R object and return a Python str

    def run(self, pdb_file_path: str, pdb_code: str) -> list:
        structure = self.pdbe.load_structure(dir=pdb_file_path, file_name=pdb_code)
        range = robjects.IntVector([1, 2, 3, 4])
        output = self.pdbe.encode_block(chain="A", range=range, structure=structure)
        print("Output is {output}".format(output=output))
        range = robjects.IntVector([5, 6, 7, 8])
        output = self.pdbe.encode_block(chain="A", range=range, structure=structure)
        print("Output is {output}".format(output=output))


if __name__ == "__main__":
    p = PDBEncode(scrap_folder="procedures")
    p.run(pdb_file_path="procedures", pdb_code="2p5u_1.pdb")
