from parsers.icsv import ICSV
import csv
import logging
import gzip


class GnomadFilter(ICSV):
    """Represents a Clinvar importer"""

    def __init__(
        self,
        gnomad_source_filename: str,
        gnomad_filtered_filename: str,
        compressed: bool,
    ):
        super().__init__(
            filename=gnomad_source_filename,
            file_type="gnomad_filter",
            max_items=-1,
            delimiter="\t",
            header_filter_prefix="##",
            compressed=compressed,
        )

        self.num_entries = 0
        self.num_snvs = 0
        self.compressed = compressed
        logging.getLogger().setLevel(logging.INFO)
        self.gnomad_filtered_filename = gnomad_filtered_filename

    def set_output_files(self):
        if self.compressed:
            self.gnomad_filter_file = gzip.open(
                filename=self.gnomad_filtered_filename, mode="wt"
            )
        else:
            self.gnomad_filter_file = open(file=self.gnomad_filtered_filename, mode="w")
        self.gnomad_filter_csv = csv.writer(
            self.gnomad_filter_file, delimiter="\t", quotechar="'"
        )

    def close(self):
        print(
            "Total entries: {num_entries}, SNVs: {num_snvs}".format(
                num_entries=self.num_entries, num_snvs=self.num_snvs
            )
        )
        self.gnomad_filter_file.close()
        super().close()

    def importfile(self):
        self.set_output_files()

        # performs standard import (calls process_item())
        super().importfile()

    def process_item(self, item: dict):
        self.num_entries += 1

        if "variant_type=snv" not in item["INFO"]:
            # skip variant if not of SNV type
            return

        info_col = item["INFO"].split(";")
        info = dict()

        # if "variant_type=snv" not in info_col:
        #     # skip variant if not of SNV type
        #     return

        for pair in info_col:
            key, value = pair.partition("=")[::2]
            info[key] = value

        # if info["variant_type"] != "snv":
        #     # skip variant if not of SNV type
        #     return

        # got record of SNV

        self.num_snvs += 1
        self.gnomad_filter_csv.writerow(item.values())
