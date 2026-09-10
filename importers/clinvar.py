from parsers.icsv import ICSV
import csv
import logging
import gzip


class ClinvarFilter(ICSV):
    """Represents a Clinvar filter"""

    def __init__(
        self,
        clinvar_source_filename: str,
        clinvar_filtered_filename: str,
        compressed: bool = False,
    ):
        super().__init__(
            filename=clinvar_source_filename,
            file_type="clinvar_filter",
            max_items=-1,
            delimiter="\t",
            header_filter_prefix="##",
            compressed=compressed,
        )
        self.num_entries = 0
        self.num_snvs = 0
        self.compressed = compressed
        logging.getLogger().setLevel(logging.INFO)
        self.clinvar_filtered_filename = clinvar_filtered_filename

    def set_output_files(self):
        if self.compressed:
            self.clinvar_filter_file = gzip.open(
                filename=self.clinvar_filtered_filename, mode="wt"
            )
        else:
            self.clinvar_filter_file = open(
                file=self.clinvar_filtered_filename, mode="w"
            )
        self.clinvar_filter_csv = csv.writer(
            self.clinvar_filter_file, delimiter="\t", quotechar="'"
        )

    def close(self):
        print(
            "Total entries: {num_entries}, SNVs: {num_snvs}".format(
                num_entries=self.num_entries, num_snvs=self.num_snvs
            )
        )
        self.clinvar_filter_file.close()
        super().close()

    def importfile(self):
        self.set_output_files()

        # performs standard import (calls process_item())
        super().importfile()

    def process_item(self, item: dict):
        self.num_entries += 1

        info_col = item["INFO"].split(";")
        info = dict()

        for pair in info_col:
            key, value = pair.partition("=")[::2]
            info[key] = value

        if info["CLNVC"] != "single_nucleotide_variant":
            # skip variant if not of SNV type
            return

        # got record of SNV
        # add chromosome prefix to the chromosome column to consistently match the chromosome format in other datasets (e.g., gnomAD)
        item["#CHROM"] = f"chr{item['#CHROM']}"

        self.num_snvs += 1
        self.clinvar_filter_csv.writerow(item.values())
