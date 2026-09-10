class FileNames:
    """Stores source and destination filenames"""

    def __init__(self, prefix: str, scrap_root: str, run_prefix=""):
        self.prefix = prefix
        source_prefix = prefix + "/source"
        output_prefix = prefix + "/output"

        self.source_prefix = source_prefix
        self.output_prefix = output_prefix

        graph_output_prefix = (
            f"{output_prefix}/graph/{run_prefix}"
            if run_prefix
            else f"{output_prefix}/graph"
        )

        self.graph_output = graph_output_prefix

        run_prefix = run_prefix + "_" if run_prefix else ""

        self.ppins_dataset_folder = f"{source_prefix}/ppins"

        self.uniprot_datafile = source_prefix + "/uniprot/2024_04/uniprot_sprot.xml"
        self.uniprot_ids_filename = (
            f"{source_prefix}/uniprot/2024_04/HUMAN_9606_idmapping_selected.tab"
        )
        self.protein2interpro_datafile = (
            f"{source_prefix}/interpro/108.0/protein2ipr.dat.gz"
        )
        self.interpro_datafile = f"{source_prefix}/interpro/108.0/interpro.xml"
        self.entrylist_datafile = f"{source_prefix}/interpro/108.0/entry.list.txt"
        self.biogrid_datafile = source_prefix + "/biogrid/BIOGRID-ALL-4.2.193.tab3.txt"
        self.clinvar_datafile = (
            source_prefix + "/clinvar/2026_06/clinvar_20260627.vcf.gz"
        )
        self.gnomad_datafile = (
            source_prefix + "/gnomad/gnomad.exomes.r2.1.1.sites.liftover_grch38.vcf.bgz"
        )
        self.sifts_datafile = (
            source_prefix + "/pdbe/sifts/uniprot_segments_observed.csv"
        )
        self.homodimers_datafile = (
            source_prefix
            + "/protein-scores-analysis/benchmark/v3/INTEGRATED_DATASET.csv"
        )

        self.clinvar_snv_filtered_filename = (
            source_prefix + "/clinvar/2026_06/clinvar_20260627.snvs.vcf.gz"
        )
        self.gnomad_snv_filtered_filename = (
            source_prefix + "/gnomad/filter_snv_gnomad.liftover_grch38.vcf.gz"
        )
        self.ferreira_etal_exome_filename = (
            source_prefix + "/ferreira_etal_exome/41586_2021_4103_MOESM4_ESM.xlsx"
        )
        self.snvs_mapped_pdb_filename = output_prefix + "/snvs_mapped_pdb.csv"
        self.snvs_mapped_pdb_dedup_filename = (
            output_prefix + "/snvs_mapped_pdb_dedup.csv"
        )
        self.vep_clinvar_filename = (
            source_prefix + "/clinvar/2026_06/vep_clinvar.vcf.gz"
        )
        self.vep_gnomad_filename = (
            source_prefix + "/gnomad/vep116_gnomad.liftover_grch38.vcf.gz"
        )
        # ---> has gaps in maps from uniprot to pdb
        # sifts_datafile = prefix + "/Databases/pdbe/sifts/pdb_chain_uniprot.csv"

        # transcriptomics (GTex)
        self.gtex_sample_attributes_filename = (
            source_prefix + "/gtex/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt"
        )
        self.gtex_transcript_tpms_filename = (
            source_prefix
            + "/gtex/GTEx_Analysis_2017-06-05_v8_RSEMv1.3.0_transcript_tpm.gct.gz"
        )

        # PDBE structures - self.pdb_database_root = source_prefix + "/pdbe/structures/divided/pdb"
        self.pdb_database_root = source_prefix + "/pdbe/biounit/pdb/divided"
        self.alphafold_database_root = (
            source_prefix + "/alphafold/UP000005640_9606_HUMAN_v4"
        )
        self.homodimers_pdb_database_root = (
            source_prefix
            + "/protein-scores-analysis/benchmark/v3/pdb/01_PDB_legacy_3dcomplex_expanded"
        )

        self.scrap_folder = f"{scrap_root}"
        self.vep_cache_dir = source_prefix + "/vep_data_116"

        self.node_protein_filename = (
            f"{graph_output_prefix}/{run_prefix}node_protein.csv"
        )
        self.node_organism_filename = (
            f"{graph_output_prefix}/{run_prefix}node_organism.csv"
        )
        self.node_interpro_entries_filename = (
            f"{graph_output_prefix}/{run_prefix}node_interpro_entries.csv"
        )
        self.rel_interpro_filename = (
            f"{graph_output_prefix}/{run_prefix}rel_interpro.csv"
        )
        self.rel_protein_interpro_entries_filename = (
            f"{graph_output_prefix}/{run_prefix}rel_protein_interpro_entries.csv"
        )
        self.node_gene_filename = f"{graph_output_prefix}/{run_prefix}node_gene.csv"
        self.node_transcript_filename = (
            f"{graph_output_prefix}/{run_prefix}node_transcript.csv"
        )
        self.node_pdb_filename = f"{graph_output_prefix}/{run_prefix}node_pdb.csv"
        self.node_model_filename = f"{graph_output_prefix}/{run_prefix}node_model.csv"
        self.node_snv_filename = f"{graph_output_prefix}/{run_prefix}node_snv.csv"
        self.node_chain_filename = f"{graph_output_prefix}/{run_prefix}node_chain.csv"
        self.node_residue_filename = (
            f"{graph_output_prefix}/{run_prefix}node_residue.csv"
        )

        self.rel_sampled_from_filename = (
            f"{graph_output_prefix}/{run_prefix}rel_sampled_from.csv"
        )
        self.rel_interacts_filename = (
            f"{graph_output_prefix}/{run_prefix}rel_interacts.csv"
        )
        self.rel_gene_transcript_filename = (
            f"{graph_output_prefix}/{run_prefix}rel_gene_transcript.csv"
        )
        self.rel_transcript_protein_filename = (
            f"{graph_output_prefix}/{run_prefix}rel_transcript_protein.csv"
        )
        self.rel_protein_pdb_filename = (
            f"{graph_output_prefix}/{run_prefix}rel_protein_pdb.csv"
        )
        self.rel_protein_model_filename = (
            f"{graph_output_prefix}/{run_prefix}rel_protein_model.csv"
        )
        self.rel_protein_snv_filename = (
            f"{graph_output_prefix}/{run_prefix}rel_protein_snv.csv"
        )
        self.rel_pdb_chain_filename = (
            f"{graph_output_prefix}/{run_prefix}rel_pdb_chain.csv"
        )
        self.rel_model_chain_filename = (
            f"{graph_output_prefix}/{run_prefix}rel_model_chain.csv"
        )
        self.rel_pdb_residue_filename = (
            f"{graph_output_prefix}/{run_prefix}rel_pdb_residue.csv"
        )
        self.rel_model_residue_filename = (
            f"{graph_output_prefix}/{run_prefix}rel_model_residue.csv"
        )
        self.rel_chain_residue_filename = (
            f"{graph_output_prefix}/{run_prefix}rel_chain_residue.csv"
        )
        self.rel_snv_residue_filename = (
            f"{graph_output_prefix}/{run_prefix}rel_snv_residue.csv"
        )
        self.rel_residue_residue_filename = (
            f"{graph_output_prefix}/{run_prefix}rel_residue_residue.csv"
        )

        self.pdb_issues_filename = f"{self.scrap_folder}/pdb_issues.csv"
        self.model_issues_filename = f"{self.scrap_folder}/model_issues.csv"
        self.physio_filename = output_prefix + "/physio.json"
        self.nonphysio_filename = output_prefix + "/nonphysio.json"

        self.pdb_residues_progress_filename = (
            f"{self.scrap_folder}/pdb_residues_progress.csv"
        )
        self.model_residues_progress_filename = (
            f"{self.scrap_folder}/model_residues_progress.csv"
        )
