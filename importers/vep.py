import parsers.vep
import time
import datetime
import csv


class VEPImporter:
    """Represents a VEP importer"""

    def __init__(self, filename):
        self.data = []
        self.parser = parsers.vep.VEPParser(filename)
        self.num_sel = 0
        self.num_read = 0
        self.num_imported = 0
        self.report_line = (
            "VEP items imported: {num_items:7d}   |   "
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

        gene_transcript_keys = set()
        gene_keys = set()
        transcript_keys = set()
        transcript_protein_keys = set()
        protein_snv_keys = set()

        # removed isoform since only the MANE_SELECT is taken
        # isoform_keys = set()

        node_gene_file = open("output/node_gene.csv", "w")
        node_gene_csv = csv.writer(node_gene_file, delimiter=";", quotechar="'")
        node_gene_csv.writerow(["id:ID(gene-id)", "symbol", "symbol_source"])

        node_transcript_file = open("output/node_transcript.csv", "w")
        node_transcript_csv = csv.writer(
            node_transcript_file, delimiter=";", quotechar="'"
        )
        node_transcript_csv.writerow(
            [
                "id:ID(transcript-id)",
            ]
        )

        rel_gene_transcript_file = open("output/rel_gene_transcript.csv", "w")
        rel_gene_transcript_csv = csv.writer(
            rel_gene_transcript_file, delimiter=";", quotechar="'"
        )
        rel_gene_transcript_csv.writerow(
            [":START_ID(gene-id)", "ref_seq", ":END_ID(transcript-id)", ":TYPE"]
        )

        node_snv_file = open("output/node_snv.csv", "w")
        node_snv_csv = csv.writer(node_snv_file, delimiter=";", quotechar="'")
        node_snv_csv.writerow(
            [
                "source_id:ID(snv-id)",
                "source",
                "allele",
                "location",
                "cDNA_position",
                "CDS_position",
                "protein_position",
                "amino_acids",
                "codons",
                "existing_variation",
                "impact",
                "sift",
                "polyphen",
                "freq",
                "clinvar_sig",
                "somatic",
                "pheno",
            ]
        )

        rel_protein_snv_file = open("output/rel_protein_snv.csv", "w")
        rel_protein_snv_csv = csv.writer(
            rel_protein_snv_file, delimiter=";", quotechar="'"
        )
        rel_protein_snv_csv.writerow(
            [":START_ID(protein-id)", ":END_ID(snv-id)", ":TYPE"]
        )

        # node_isoform_file = open("node_isoform.csv", "w")
        # node_isoform_csv = csv.writer(node_isoform_file, delimiter=";", quotechar="'")
        # node_isoform_csv.writerow(
        #     [
        #         "accession:ID(isoform-id)",
        #     ]
        # )

        rel_transcript_protein_file = open("output/rel_transcript_protein.csv", "w")
        rel_transcript_protein_csv = csv.writer(
            rel_transcript_protein_file, delimiter=";", quotechar="'"
        )
        rel_transcript_protein_csv.writerow(
            [":START_ID(transcript-id)", "version", ":END_ID(protein-id)", ":TYPE"]
        )

        # rel_transcript_isoform_file = open("rel_transcript_isoform.csv", "w")
        # rel_transcript_isoform_csv = csv.writer(
        #     rel_transcript_isoform_file, delimiter=";", quotechar="'"
        # )
        # rel_transcript_isoform_csv.writerow(
        #     [":START_ID(transcript-id)", ":END_ID(isoform-id)", ":TYPE"]
        # )

        while True:
            item = self.parser.get_next_item()

            self.num_read += 1

            if item is None:
                print("\nNo more items")
                node_gene_file.close()
                node_transcript_file.close()
                rel_gene_transcript_file.close()
                node_snv_file.close()
                rel_protein_snv_file.close()
                # node_isoform_file.close()
                rel_transcript_protein_file.close()
                # rel_transcript_isoform_file.close()
                break

            extra_col = item["Extra"].split(";")
            extra = dict()

            for pair in extra_col:
                key, value = pair.partition("=")[::2]
                extra[key] = value

            # basic filtering
            if "SWISSPROT" not in extra:
                continue

            if "MANE_SELECT" not in extra:
                continue

            if item["Feature_type"] != "Transcript":
                continue

            # -> missense_variant can come as multivalued consequence (e.g. "missense_variant,splice_region_variant")
            if "missense_variant" not in item["Consequence"]:
                continue

            # GENE node
            gene = dict()
            gene["id"] = item["Gene"]
            gene["symbol"] = ""
            if "SYMBOL" in extra:
                gene["symbol"] = extra["SYMBOL"]
            gene["symbol_source"] = ""
            if "SYMBOL_SOURCE" in extra:
                gene["symbol_source"] = extra["SYMBOL_SOURCE"]
            if gene["id"] not in gene_keys:
                node_gene_csv.writerow(gene.values())
                gene_keys.add(gene["id"])

            # TRANSCRIPT node
            transcript = dict()
            transcript["id"] = item["Feature"]
            if transcript["id"] not in transcript_keys:
                node_transcript_csv.writerow(transcript.values())
                transcript_keys.add(transcript["id"])

            # GENE->TRANSCRIPT relation
            gene_transcript_key = gene["id"] + extra["MANE_SELECT"] + transcript["id"]
            if gene_transcript_key not in gene_transcript_keys:
                rel_gene_transcript_csv.writerow(
                    [gene["id"], extra["MANE_SELECT"], transcript["id"], "reads"]
                )
                gene_transcript_keys.add(gene_transcript_key)

            # SNV node
            swissprot = extra["SWISSPROT"].split(".")
            swissprot_id = swissprot[0]
            swissprot_version = swissprot[1]

            snv = dict()
            snv["source_id"] = item["#Uploaded_variation"]
            snv["source"] = "clinvar"
            snv["allele"] = item["Allele"]
            snv["location"] = item["Location"]
            snv["cDNA_position"] = item["cDNA_position"]
            snv["CDS_position"] = item["CDS_position"]
            snv["protein_position"] = item["Protein_position"]
            snv["amino_acids"] = item["Amino_acids"]
            snv["codons"] = item["Codons"]
            snv["existing_variation"] = item["Existing_variation"]
            snv["impact"] = "-"
            if "IMPACT" in extra:
                snv["impact"] = extra["IMPACT"]
            snv["sift"] = "-"
            if "SIFT" in extra:
                snv["sift"] = extra["SIFT"]
            snv["polyphen"] = "-"
            if "PolyPhen" in extra:
                snv["polyphen"] = extra["PolyPhen"]
            snv["freq"] = "-"
            if "AF" in extra:
                snv["freq"] = extra["AF"]
            snv["clinvar_sig"] = "-"
            if "CLIN_SIG" in extra:
                snv["clinvar_sig"] = extra["CLIN_SIG"]
            snv["somatic"] = "-"
            if "SOMATIC" in extra:
                snv["somatic"] = extra["SOMATIC"]
            snv["pheno"] = "-"
            if "PHENO" in extra:
                snv["pheno"] = extra["PHENO"]

            node_snv_csv.writerow(snv.values())

            # PROTEIN->SNV relation
            protein_snv_key = swissprot_id + snv["source_id"]
            if protein_snv_key not in protein_snv_keys:
                rel_protein_snv_csv.writerow(
                    [swissprot_id, snv["source_id"], "contains"]
                )
                protein_snv_keys.add(protein_snv_key)

            # ISOFORM node - removed since only the MANE_SELECT is taken
            # isoform = dict()
            # isoform["id"] = swissprot_id
            # if "UNIPROT_ISOFORM" in extra:
            #     isoform["id"] = extra["UNIPROT_ISOFORM"]
            # if isoform["id"] not in isoform_keys:
            #     node_isoform_csv.writerow(isoform.values())
            #     isoform_keys.add(isoform["id"])

            # TRANSCRIPT->PROTEIN relation
            transcript_protein_key = transcript["id"] + swissprot_version + swissprot_id
            if transcript_protein_key not in transcript_protein_keys:
                rel_transcript_protein_csv.writerow(
                    [transcript["id"], swissprot_version, swissprot_id, "generates"]
                )
                transcript_protein_keys.add(transcript_protein_key)

            # TRANSCRIPT->ISOFORM relation
            # rel_transcript_isoform_csv.writerow(
            #     [transcript["id"], isoform["id"], "maps"]
            # )

            self.num_imported += 1
            if self.num_read % 1000 == 0:
                stop_counter = time.perf_counter()
                speed = 1000 / (stop_counter - start_counter)
                eta = str(datetime.timedelta(seconds=round(11413278 / speed)))
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
