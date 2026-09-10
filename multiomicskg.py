import logging
import gc
from os import makedirs

from importers.uniprot import UniprotImporter
from importers.ppins.ppins import PPINs

from importers.interpro import InterproImporter
from importers.clinvar import ClinvarFilter
from importers.gnomad import GnomadFilter

from importers.sifts import SIFTSImporter

from importers.snv_residues import SNVResidues
from importers.vep_snv import VepSnv
from importers.gtex_transcript import GTexTranscript

import argparse

from procedures.vep_caller import VepCaller

import filenames

if __name__ == "__main__":
    logging.getLogger().setLevel(logging.INFO)

    parser = argparse.ArgumentParser(
        description="Multiomics database integration and generation tool",
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
        "-v",
        required=True,
        type=str,
        default=False,
        dest="vep_path",
        help="path of vep tool, e.g.: ./deps/ensembl-vep/vep",
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

    fn = filenames.FileNames(args.datastore_prefix, "", args.run_prefix)
    vep_path = args.vep_path

    makedirs(fn.graph_output, exist_ok=True)

    # reads: uniprot_sprot.xml
    # writes: node_protein.csv, node_organism.csv, rel_sampled_from.csv, node_model.csv, rel_protein_model.csv
    uniprotimp = UniprotImporter(
        fn.uniprot_datafile,
        fn.node_protein_filename,
        fn.node_organism_filename,
        fn.rel_sampled_from_filename,
        fn.node_model_filename,
        fn.rel_protein_model_filename,
        human_only=True,
    )
    uniprotimp.importdb()
    uniprotimp.close()
    del uniprotimp
    gc.collect()

    # reads: protein2ipr.dat.gz, interpro.xml, entry.list.txt, node_protein.csv
    # writes: node_interpro_entries.csv, rel_protein_interpro_entries.csv, rel_interpro.csv
    interproimp = InterproImporter(
        protein2interpro_filename=fn.protein2interpro_datafile,
        interpro_filename=fn.interpro_datafile,
        entrylist_filename=fn.entrylist_datafile,
        node_interpro_entries_filename=fn.node_interpro_entries_filename,
        rel_interpro_filename=fn.rel_interpro_filename,
        rel_protein_interpro_entries_filename=fn.rel_protein_interpro_entries_filename,
        node_protein_filename=fn.node_protein_filename,
    )
    interproimp.importdb()
    del interproimp
    gc.collect()

    # reads: node_protein.csv, ppins dataset files, HUMAN_9606_idmapping_selected.tab
    # writes: rel_interacts.csv
    ppins = PPINs(
        node_protein_filename=fn.node_protein_filename,
        ppins_dataset_folder=fn.ppins_dataset_folder,
        rel_interacts_filename=fn.rel_interacts_filename,
        uniprot_ids_filename=fn.uniprot_ids_filename,
    )
    ppins.create_unified_ppin()
    del ppins
    gc.collect()

    # reads: clinvar_20260627.vcf.gz
    # writes: clinvar_20260627.snvs.vcf.gz
    clinvarfilter = ClinvarFilter(
        clinvar_source_filename=fn.clinvar_datafile,
        clinvar_filtered_filename=fn.clinvar_snv_filtered_filename,
        compressed=True,
    )
    clinvarfilter.importfile()
    clinvarfilter.close()
    del clinvarfilter
    gc.collect()

    # reads: gnomad.exomes.r2.1.1.sites.liftover_grch38.vcf.bgz
    # writes: filter_snv_gnomad.liftover_grch38.vcf.gz
    gnomadfilter = GnomadFilter(
        gnomad_source_filename=fn.gnomad_datafile,
        gnomad_filtered_filename=fn.gnomad_snv_filtered_filename,
        compressed=True,
    )
    gnomadfilter.importfile()
    gnomadfilter.close()
    del gnomadfilter
    gc.collect()

    # reads: clinvar_20260627.snvs.vcf.gz
    # writes: vep_clinvar.vcf.gz
    vep_clinvar = VepCaller(
        source_filename=fn.clinvar_snv_filtered_filename,
        output_filename=fn.vep_clinvar_filename,
        vep_path=vep_path,
        vep_cache_dir=fn.vep_cache_dir,
        clinvar_source_filename=fn.clinvar_datafile,
    )
    if vep_clinvar.run() is False:
        exit(1)

    # reads: filter_snv_gnomad.liftover_grch38.vcf.gz
    # writes: vep_gnomad.liftover_grch38.vcf.gz
    vep_gnomad = VepCaller(
        source_filename=fn.gnomad_snv_filtered_filename,
        output_filename=fn.vep_gnomad_filename,
        vep_path=vep_path,
        vep_cache_dir=fn.vep_cache_dir,
    )
    if vep_gnomad.run() is False:
        exit(1)

    # reads: vep_clinvar.vcf.gz, vep116_gnomad.liftover_grch38.vcf.gz, 41586_2021_4103_MOESM4_ESM.xlsx
    # writes: node_gene.csv, node_transcript.csv, rel_gene_transcript.csv, node_snv.csv, rel_protein_snv.csv, rel_transcript_protein.csv
    # vepimp = VEPImporter(vep_datafile)
    # vepimp.importdb()
    # vepimp.close()
    vep_snv = VepSnv(
        vep_clinvar_filename=fn.vep_clinvar_filename,
        vep_gnomad_filename=fn.vep_gnomad_filename,
        ferreira_etal_exome_filename=fn.ferreira_etal_exome_filename,
        node_gene_filename=fn.node_gene_filename,
        node_transcript_filename=fn.node_transcript_filename,
        rel_gene_transcript_filename=fn.rel_gene_transcript_filename,
        node_snv_filename=fn.node_snv_filename,
        rel_protein_snv_filename=fn.rel_protein_snv_filename,
        rel_transcript_protein_filename=fn.rel_transcript_protein_filename,
    )
    vep_snv.importfile()
    vep_snv.close()
    del vep_snv
    gc.collect()

    # # reads: GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt, GTEx_Analysis_2017-06-05_v8_RSEMv1.3.0_transcript_tpm.gct.gz
    # # writes: node_transcript.csv
    gtex_transcripts = GTexTranscript(
        gtex_sample_attributes_filename=fn.gtex_sample_attributes_filename,
        gtex_transcript_tpms_filename=fn.gtex_transcript_tpms_filename,
        node_transcript_filename=fn.node_transcript_filename,
    )
    gtex_transcripts.importfile()

    # reads: uniprot_segments_observed.csv (SIFTS dataset), node_protein.csv (optional)
    # writes: node_pdb.csv rel_protein_pdb.csv
    siftsimp = SIFTSImporter(
        filename=fn.sifts_datafile,
        node_protein_filename=fn.node_protein_filename,
        node_pdb_filename=fn.node_pdb_filename,
        rel_protein_pdb_filename=fn.rel_protein_pdb_filename,
        filter_by_node_protein=True,
    )
    siftsimp.importfile()
    siftsimp.close()
    del siftsimp
    gc.collect()

    # only for entries that map to a protein with a link to a pdb file
    # reads: node_snv.csv, rel_protein_snv.csv, rel_protein_pdb.csv, rel_protein_model.csv
    # writes: rel_snv_residue.csv
    snv_residues = SNVResidues(
        node_snv_filename=fn.node_snv_filename,
        rel_protein_snv_filename=fn.rel_protein_snv_filename,
        rel_protein_pdb_filename=fn.rel_protein_pdb_filename,
        rel_protein_model_filename=fn.rel_protein_model_filename,
        rel_snv_residue_filename=fn.rel_snv_residue_filename,
    )
    snv_residues.importfile()
    snv_residues.close()
    del snv_residues
    gc.collect()
