import subprocess
import logging
import multiprocessing


class VepCaller:
    """Executes a VEP session from a list of variants in VCF"""

    def __init__(
        self,
        source_filename: str,
        output_filename: str,
        vep_path: str,
        vep_cache_dir: str,
        clinvar_source_filename: str = "",
    ):
        self.source_filename = source_filename
        self.output_filename = output_filename
        self.vep_path = vep_path
        self.vep_cache_dir = vep_cache_dir
        self.clinvar_source_filename = clinvar_source_filename
        logging.getLogger().setLevel(logging.INFO)

    def run(self) -> bool:
        try:
            num_cpus = multiprocessing.cpu_count() - 1
            logging.info("Running VEP on " + self.source_filename)
            args = [
                "./deps/ensembl-vep/vep",
                "--af",
                "--af_gnomad",
                "-fork",
                str(num_cpus),
                "--appris",
                "--biotype",
                "--buffer_size",
                "5000",
                "--check_existing",
                "--distance",
                "5000",
                "--mane",
                "--polyphen",
                "b",
                "--protein",
                "--pubmed",
                "--regulatory",
                "--sift",
                "b",
                "--species",
                "homo_sapiens",
                "--symbol",
                "--transcript_version",
                "--tsl",
                "--uniprot",
                "--cache",
                "--input_file",
                self.source_filename,
                "--output_file",
                self.output_filename,
                "--offline",
                "--force_overwrite",
                "--verbose",
                "--dir",
                self.vep_cache_dir,
                "--compress_output",
                "gzip",
            ]
            if self.clinvar_source_filename:
                args.append("--clinvar_somatic_classification")
                args.append("--custom")
                args.append(f"file={self.clinvar_source_filename},short_name=ClinVar,format=vcf,type=exact,coords=0,fields=CLNSIG%CLNREVSTAT%CLNDN")

            result = subprocess.run(
                args=args, check=True, text=True, capture_output=False
            )
            if result.returncode == 0:
                logging.info("Finished VepCaller successfully")
                return True
            else:
                logging.error("Error whilst running VepCaller")
                return False
        except Exception as e:
            logging.error("Error running VepCaller: " + str(e))
            return False
