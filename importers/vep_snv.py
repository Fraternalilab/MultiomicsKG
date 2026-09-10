from parsers.icsv import ICSV
import logging
import csv
import re
import pandas as pd


class VepSnv:
    def __init__(
        self,
        vep_clinvar_filename: str,
        vep_gnomad_filename: str,
        ferreira_etal_exome_filename: str,
        node_gene_filename: str,
        node_transcript_filename: str,
        rel_gene_transcript_filename: str,
        node_snv_filename: str,
        rel_protein_snv_filename: str,
        rel_transcript_protein_filename: str,
    ):
        self.vep_clinvar = ICSV(
            filename=vep_clinvar_filename,
            file_type="vep_clinvar",
            delimiter="\t",
            header_filter_prefix="##",
            max_items=-1,
            compressed=True,
        )
        self.vep_gnomad = ICSV(
            filename=vep_gnomad_filename,
            file_type="vep_gnomad",
            delimiter="\t",
            header_filter_prefix="##",
            max_items=-1,
            compressed=True,
        )
        self.ferreira_etal_exome_df = pd.read_excel(
            open(ferreira_etal_exome_filename, "rb"), "ST6", header=1
        )
        self.node_gene_filename = node_gene_filename
        self.node_transcript_filename = node_transcript_filename
        self.rel_gene_transcript_filename = rel_gene_transcript_filename
        self.node_snv_filename = node_snv_filename
        self.rel_protein_snv_filename = rel_protein_snv_filename
        self.rel_transcript_protein_filename = rel_transcript_protein_filename
        self.num_items = 0
        self.polyphen_sifts_pattern = re.compile("(.+)\((.+)\)")
        logging.getLogger().setLevel(logging.INFO)

    def close(self):
        print("\nVEP items imported: {num_items}".format(num_items=self.num_items))
        self.vep_clinvar.close()
        self.vep_gnomad.close()
        self.node_gene_file.close()
        self.node_transcript_file.close()
        self.rel_gene_transcript_file.close()
        self.node_snv_file.close()
        self.rel_protein_snv_file.close()
        self.rel_transcript_protein_file.close()

    def set_output_files(self):
        self.node_gene_file = open(self.node_gene_filename, "w")
        self.node_gene_csv = csv.writer(
            self.node_gene_file,
            delimiter=";",
            quotechar="'",
            escapechar="\\",
            quoting=csv.QUOTE_ALL,
        )
        self.node_gene_csv.writerow(
            [
                "id:ID(gene-id)",
                "name:string",
                "symbol_source:string",
                "trait:string",
                "trait_description:string",
                "ukb_detailed_trait_name:string",
            ]
        )

        self.node_transcript_file = open(self.node_transcript_filename, "w")
        self.node_transcript_csv = csv.writer(
            self.node_transcript_file, delimiter=";", quotechar="'"
        )
        self.node_transcript_csv.writerow(
            [
                "id:ID(transcript-id)",
            ]
        )

        self.rel_gene_transcript_file = open(self.rel_gene_transcript_filename, "w")
        self.rel_gene_transcript_csv = csv.writer(
            self.rel_gene_transcript_file, delimiter=";", quotechar="'"
        )
        self.rel_gene_transcript_csv.writerow(
            [":START_ID(gene-id)", "ref_seq", ":END_ID(transcript-id)", ":TYPE"]
        )

        self.node_snv_file = open(self.node_snv_filename, "w")
        self.node_snv_csv = csv.writer(self.node_snv_file, delimiter=";", quotechar="'")
        self.node_snv_csv.writerow(
            [
                "id:ID(snv-id)",
                "source_id:string",
                "source:string",
                "allele:string",
                "location:string",
                "cDNA_position:long",
                "CDS_position:long",
                "protein_position:long",
                "amino_acids:string",
                "codons:string",
                "lead_marker:string",
                "existing_variation:string",
                "impact:string",
                "mane_select:string",
                "effect_direction:string",
                "effect_p-value:double",
                "effect_allele_frequency:double",
                "effect_(95%_ci):string",
                "sift:double",
                "sift_significance:string",
                "polyphen:double",
                "polyphen_significance:string",
                "af_1kg:double",
                "gnomad_af:double",
                "gnomad_afr_af:double",
                "gnomad_amr_af:double",
                "gnomad_asj_af:double",
                "gnomad_eas_af:double",
                "gnomad_fin_af:double",
                "gnomad_mid_af:double",
                "gnomad_nfe_af:double",
                "gnomad_rem_af:double",
                "gnomad_sas_af:double",
                "clinvar_all_significances:string",
                "clinvar_id:string",
                "clinvar_disease_name:string",
                "clinvar_review_status:string",
                "clinvar_significance:string",
                "somatic:string",
                "pheno:string",
            ]
        )

        self.rel_protein_snv_file = open(self.rel_protein_snv_filename, "w")
        self.rel_protein_snv_csv = csv.writer(
            self.rel_protein_snv_file, delimiter=";", quotechar="'"
        )
        self.rel_protein_snv_csv.writerow(
            [":START_ID(protein-id)", ":END_ID(snv-id)", ":TYPE"]
        )

        # node_isoform_file = open("node_isoform.csv", "w")
        # node_isoform_csv = csv.writer(node_isoform_file, delimiter=";", quotechar="'")
        # node_isoform_csv.writerow(
        #     [
        #         "accession:ID(isoform-id)",
        #     ]
        # )

        self.rel_transcript_protein_file = open(
            self.rel_transcript_protein_filename, "w"
        )
        self.rel_transcript_protein_csv = csv.writer(
            self.rel_transcript_protein_file, delimiter=";", quotechar="'"
        )
        self.rel_transcript_protein_csv.writerow(
            [":START_ID(transcript-id)", "version", ":END_ID(protein-id)", ":TYPE"]
        )

        # rel_transcript_isoform_file = open("rel_transcript_isoform.csv", "w")
        # rel_transcript_isoform_csv = csv.writer(
        #     rel_transcript_isoform_file, delimiter=";", quotechar="'"
        # )
        # rel_transcript_isoform_csv.writerow(
        #     [":START_ID(transcript-id)", ":END_ID(isoform-id)", ":TYPE"]
        # )

    def importfile(self):
        self.set_output_files()

        self.gene_transcript_keys = set()
        self.gene_keys = set()
        self.transcript_keys = set()
        self.transcript_protein_keys = set()
        self.protein_snv_keys = set()

        # removed isoform since only the MANE_SELECT is taken
        # isoform_keys = set()

        # load clinvar into memory
        self.clinvar_db = self.load_clinvar()

        while True:
            # process gnomad item
            item = self.vep_gnomad.get_next_item()

            if item is None:
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

            item["Extra"] = extra
            item["source_db"] = "gnomad"

            key = "{location}_{allele}_{feature}_{aminoacids}_{codons}".format(
                location=item["Location"],
                allele=item["Allele"],
                feature=item["Feature"],
                aminoacids=item["Amino_acids"],
                codons=item["Codons"],
            )

            if key in self.clinvar_db:
                # merge with clinvar item
                self.merge_item(key, item)
            else:
                self.write_item(key, item)

        # dump remaining clinvar items including merged ones
        for key, item in self.clinvar_db.items():
            self.write_item(key, item)
        print("Finished")

    def load_clinvar(self) -> dict:
        clinvar_db = dict()

        while True:
            item = self.vep_clinvar.get_next_item()

            if item is None:
                return clinvar_db

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

            item["Extra"] = extra
            item["source_db"] = "clinvar"
            item["merged"] = False

            key = "{location}_{allele}_{feature}_{aminoacids}_{codons}".format(
                location=item["Location"],
                allele=item["Allele"],
                feature=item["Feature"],
                aminoacids=item["Amino_acids"],
                codons=item["Codons"],
            )
            clinvar_db[key] = item

    def merge_item(self, key: str, gnomad_item: dict):
        if self.clinvar_db[key]["merged"]:
            print("Item already merged: " + key)
            return

        # add missing fields to gnomad_item so the comparison works
        if "ClinVar" in self.clinvar_db[key]["Extra"]:
            gnomad_item["Extra"]["ClinVar"] = self.clinvar_db[key]["Extra"]["ClinVar"]
        if "ClinVar_CLNDN" in self.clinvar_db[key]["Extra"]:
            gnomad_item["Extra"]["ClinVar_CLNDN"] = self.clinvar_db[key]["Extra"][
                "ClinVar_CLNDN"
            ]
        if "ClinVar_CLNREVSTAT" in self.clinvar_db[key]["Extra"]:
            gnomad_item["Extra"]["ClinVar_CLNREVSTAT"] = self.clinvar_db[key]["Extra"][
                "ClinVar_CLNREVSTAT"
            ]
        if "ClinVar_CLNSIG" in self.clinvar_db[key]["Extra"]:
            gnomad_item["Extra"]["ClinVar_CLNSIG"] = self.clinvar_db[key]["Extra"][
                "ClinVar_CLNSIG"
            ]
        if "ClinVar_FILTER" in self.clinvar_db[key]["Extra"]:
            gnomad_item["Extra"]["ClinVar_FILTER"] = self.clinvar_db[key]["Extra"][
                "ClinVar_FILTER"
            ]

        added, removed, modified, same = self.compare_dict(
            self.clinvar_db[key], gnomad_item
        )
        if len(modified) == 2:
            if "#Uploaded_variation" in modified and "source_db" in modified:
                self.clinvar_db[key]["#Uploaded_variation"] = (
                    self.clinvar_db[key]["#Uploaded_variation"]
                    + ","
                    + gnomad_item["#Uploaded_variation"]
                )
                self.clinvar_db[key]["merged"] = True
                self.clinvar_db[key]["source_db"] = "clinvar,gnomad"
                return

        print("merge here!!!")
        print(key)
        print(self.clinvar_db[key])
        print(gnomad_item)
        return

    def compare_dict(self, d1: dict, d2: dict):
        d1_keys = set(d1.keys())
        d2_keys = set(d2.keys())
        shared_keys = d1_keys.intersection(d2_keys)
        added = d1_keys - d2_keys
        removed = d2_keys - d1_keys
        modified = {o: (d1[o], d2[o]) for o in shared_keys if d1[o] != d2[o]}
        same = set(o for o in shared_keys if d1[o] == d2[o])
        return added, removed, modified, same

    def write_item(self, key: str, item: dict):
        # GENE node
        gene = dict()
        gene["id"] = item["Gene"]
        gene["name"] = item["Extra"].get("SYMBOL", "")
        if gene["id"] not in self.gene_keys:
            gene["symbol_source"] = item["Extra"].get("SYMBOL_SOURCE", "")
            gene["trait"] = ""
            gene["trait_description"] = ""
            gene["ukb_detailed_trait_name"] = ""
            # find gene in ferreira_etal
            gene_info_s = self.ferreira_etal_exome_df.loc[
                self.ferreira_etal_exome_df["Gene"] == gene["name"]
            ]
            if gene_info_s.empty is not True:
                if gene_info_s.shape[0] != 1:
                    raise Exception("Unexpected shape of gene data frame")
                gene_info_s.index = [0]
                gene["trait"] = gene_info_s.loc[0].at["Trait"]
                gene["trait_description"] = gene_info_s.loc[0].at["Trait description"]
                gene["ukb_detailed_trait_name"] = gene_info_s.loc[0].at[
                    "UKB detailed trait name"
                ]
            self.node_gene_csv.writerow(gene.values())
            self.gene_keys.add(gene["id"])

        # TRANSCRIPT node
        transcript = dict()
        transcript["id"] = item["Feature"]
        if transcript["id"] not in self.transcript_keys:
            self.node_transcript_csv.writerow(transcript.values())
            self.transcript_keys.add(transcript["id"])

        # GENE->TRANSCRIPT relation
        gene_transcript_key = (
            gene["id"] + item["Extra"]["MANE_SELECT"] + transcript["id"]
        )
        if gene_transcript_key not in self.gene_transcript_keys:
            self.rel_gene_transcript_csv.writerow(
                [gene["id"], item["Extra"]["MANE_SELECT"], transcript["id"], "reads"]
            )
            self.gene_transcript_keys.add(gene_transcript_key)

        # SNV node
        swissprot = item["Extra"]["SWISSPROT"].split(".")
        swissprot_id = swissprot[0]
        swissprot_version = swissprot[1]

        # retrieve again variant data from ferreira_etal
        ferreira_etal = dict()
        gene_info_s = self.ferreira_etal_exome_df.loc[
            self.ferreira_etal_exome_df["Gene"] == gene["name"]
        ]
        if gene_info_s.empty is not True:
            if gene_info_s.shape[0] != 1:
                raise Exception("Unexpected shape of gene data frame")
            gene_info_s.index = [0]
            if gene_info_s.loc[0].at["Marker type"] == "DelMissense":
                location = f"chr{gene_info_s.loc[0].at['Chr']:.0f}:{gene_info_s.loc[0].at['Position']:.0f}"
                if location == item["Location"]:
                    ferreira_etal["lead_marker"] = gene_info_s.loc[0].at["Lead marker"]
                    ferreira_etal["effect_direction"] = gene_info_s.loc[0].at[
                        "Effect direction"
                    ]
                    ferreira_etal["effect_p-value"] = gene_info_s.loc[0].at["P-value"]
                    ferreira_etal["effect_allele_frequency"] = gene_info_s.loc[0].at[
                        "Effect allele frequency"
                    ]
                    ferreira_etal["effect_(95%_ci)"] = gene_info_s.loc[0].at[
                        "Effect (95% CI)"
                    ]

        # parse polyphen and sift values
        polyphen_val = item["Extra"].get("PolyPhen", "-")
        sift_val = item["Extra"].get("SIFT", "-")
        polyphen_significance = ""
        sift_significance = ""
        polyphen_factor = ""
        sift_factor = ""

        match = self.polyphen_sifts_pattern.match(polyphen_val)
        if match is not None:
            polyphen_significance = match.group(1)
            polyphen_factor = float(match.group(2))
        match = self.polyphen_sifts_pattern.match(sift_val)
        if match is not None:
            sift_significance = match.group(1)
            sift_factor = float(match.group(2))

        snv = dict()
        snv["id"] = key
        snv["source_id"] = item["#Uploaded_variation"]
        snv["source"] = item["source_db"]
        snv["allele"] = item["Allele"]
        snv["location"] = item["Location"]
        snv["cDNA_position"] = item["cDNA_position"]
        snv["CDS_position"] = item["CDS_position"]
        snv["protein_position"] = item["Protein_position"]
        snv["amino_acids"] = item["Amino_acids"]
        snv["codons"] = item["Codons"]
        snv["lead_marker"] = ferreira_etal.get("lead_marker", "")
        snv["existing_variation"] = item["Existing_variation"]
        snv["impact"] = item["Extra"].get("IMPACT", "")
        snv["mane_select"] = item["Extra"].get("MANE_SELECT", "")
        snv["effect_direction"] = ferreira_etal.get("effect_direction", "")
        snv["effect_p-value"] = ferreira_etal.get("effect_p-value", 0.0)
        snv["effect_allele_frequency"] = ferreira_etal.get(
            "effect_allele_frequency", 0.0
        )
        snv["effect_(95%_ci)"] = ferreira_etal.get("effect_(95%_ci)", "")
        snv["sift"] = sift_factor
        snv["sift_significance"] = sift_significance
        snv["polyphen"] = polyphen_factor
        snv["polyphen_significance"] = polyphen_significance
        snv["af_1kg"] = item["Extra"].get("AF", "")
        snv["gnomad_af"] = item["Extra"].get("gnomADe_AF", "")
        snv["gnomad_afr_af"] = item["Extra"].get("gnomADe_AFR_AF", "")
        snv["gnomad_amr_af"] = item["Extra"].get("gnomADe_AMR_AF", "")
        snv["gnomad_asj_af"] = item["Extra"].get("gnomADe_ASJ_AF", "")
        snv["gnomad_eas_af"] = item["Extra"].get("gnomADe_EAS_AF", "")
        snv["gnomad_fin_af"] = item["Extra"].get("gnomADe_FIN_AF", "")
        snv["gnomad_mid_af"] = item["Extra"].get("gnomADe_MID_AF", "")
        snv["gnomad_nfe_af"] = item["Extra"].get("gnomADe_NFE_AF", "")
        snv["gnomad_rem_af"] = item["Extra"].get("gnomADe_REMAINING_AF", "")
        snv["gnomad_sas_af"] = item["Extra"].get("gnomADe_SAS_AF", "")
        snv["clinvar_all_significances"] = item["Extra"].get("CLIN_SIG", "")
        snv["clinvar_id"] = item["Extra"].get("ClinVar", "")
        snv["clinvar_disease_name"] = item["Extra"].get("ClinVar_CLNDN", "")
        snv["clinvar_review_status"] = item["Extra"].get("ClinVar_CLNREVSTAT", "")
        snv["clinvar_significance"] = item["Extra"].get("ClinVar_CLNSIG", "")
        snv["somatic"] = item["Extra"].get("SOMATIC", "")
        snv["pheno"] = item["Extra"].get("PHENO", "")

        self.node_snv_csv.writerow(snv.values())

        # PROTEIN->SNV relation
        protein_snv_key = swissprot_id + snv["id"]
        if protein_snv_key not in self.protein_snv_keys:
            self.rel_protein_snv_csv.writerow([swissprot_id, snv["id"], "contains"])
            self.protein_snv_keys.add(protein_snv_key)

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
        if transcript_protein_key not in self.transcript_protein_keys:
            self.rel_transcript_protein_csv.writerow(
                [transcript["id"], swissprot_version, swissprot_id, "generates"]
            )
            self.transcript_protein_keys.add(transcript_protein_key)

        # TRANSCRIPT->ISOFORM relation
        # rel_transcript_isoform_csv.writerow(
        #     [transcript["id"], isoform["id"], "maps"]
        # )
