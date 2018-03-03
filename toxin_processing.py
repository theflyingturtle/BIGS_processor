#!/usr/bin/env python3

import pandas as pd
import xlrd
import click
import collections
import logging
import os
import shutil
import pathlib
import sys

from Bio import SeqIO
from Bio.Data.CodonTable import TranslationError


def clusters(df):
    clusters = collections.defaultdict(list)
    for column in df.columns:
        if column.startswith("mj_"):
            clusters[column.replace("mj_", "")[:3]].append(column)

    return clusters


def read_export(allele_export):
    df = pd.read_excel(allele_export)
    loci = [locus for loci in clusters(df).values() for locus in loci]
    # TO DO: script failed when I used a dataset that had no [S][I] or [S] to replace??
    # Replace "[S][I]" and "[S]" with "EOC" (end of contig) and "Q" (questionable), respectively
    df.loc[:, loci] = df.loc[:, loci].replace(["[S][I]", "[S]"], ["EOC", "Q"])
    # Drop rows where loci appear to be uncurated
    missing_data = df[df[loci].isnull().any(axis="columns")]["id"]
    if not missing_data.empty:
        logging.warn("Dropped uncurated isolate(s) with id(s): %s", ", ".join(missing_data.astype(str)))
        df.drop(missing_data.index, inplace=True)

    return df


def read_ref_seqs(ref_seqs):
    ref_lens = {}
    for record in SeqIO.parse(ref_seqs, "fasta"):
        # reference exports have _1 in their names but we just want the locus
        ref_lens[record.name.replace("_1", "")] = len(record.seq)
    
    return ref_lens


def read_isolate_alleles(allele_seqs):
    loci = {}
    for f in pathlib.Path(allele_seqs).glob("*.fas"):
        alleles = {}
        for record in SeqIO.parse(str(f), "fasta"):  # BioPython poops because it checks for type str
            alleles[int(record.name)] = record.seq
        loci[f.stem] = alleles

    return loci


def prepare_outdir(output_directory, overwrite=False):
    output_directory = pathlib.Path(output_directory)
    if output_directory.exists():
        if overwrite:
            shutil.rmtree(output_directory)
        else:
            logging.fatal("%s already exists but you didn't pass --overwrite", output_directory)
            assert False, "Cannot create output directory"

    output_directory.mkdir()
    return output_directory


@click.command()
@click.argument("allele_export", type=click.File("rb"))
@click.argument("ref_seqs", type=click.Path(exists=True, dir_okay=False))
@click.argument("allele_seqdir", type=click.Path(exists=True, file_okay=False))
@click.option("--output_directory", 
    default="outdir",
    prompt=True,
    type=click.Path(file_okay=False, dir_okay=True, readable=True, writable=True, resolve_path=True))
@click.option('--overwrite', is_flag=True)
def main(allele_export, ref_seqs, allele_seqdir, output_directory, overwrite):
    """TO DO... Short description of what the script does

    Long description of what the script does"""
    outdir = prepare_outdir(output_directory, overwrite)

    df = read_export(allele_export)

    c = clusters(df)
    logging.debug("Clusters: %s", c)

    # Summary allele presence across clusters
    cluster_summary = pd.DataFrame(index=df.index)
    for cluster, loci in c.items():
        cluster_summary[cluster] = (df[loci] != 0).sum(axis="columns") / len(loci)
    cluster_summary.to_csv(outdir / "cluster_presence_absence_summary.csv")

    # Summary allele presence within clusters
    gene_list = []
    for cluster, loci in c.items():
          gene = (df[loci] != 0).sum(axis="rows") / len(df.index)
          gene_list.append(gene)
    gene_list = pd.concat(gene_list)
    gene_prev = pd.DataFrame({'locus':gene_list.index, 'proportion_present':gene_list.values})
    gene_prev['cluster'] = gene_prev['locus'].replace('mj_', '', regex=True).str[:3]
    gene_prev = gene_prev.reindex(sorted(gene_prev.columns), axis=1)
    gene_prev.to_csv(outdir / "gene_prevalence_summary.csv")

    # To process one file
    ref_lens = read_ref_seqs(ref_seqs)
    statuses = {}
    for locus, alleles in read_isolate_alleles(allele_seqdir).items():
        status = {}
        for isolate_id, seq in alleles.items():
            # Absent == sequence starts with a gap (-)
            if seq.startswith("-"):
                status[isolate_id] = "absent"
                continue

            # Incomplete == sequence does not translate
            try:
                seq.translate(cds=True)
            except TranslationError:
                status[isolate_id] = "incomplete"
                continue

            # Fragment == sequence < 90% of reference allele length
            if len(seq) / ref_lens[locus] < 0.9:
                status[isolate_id] = "fragment"
                continue

            # Overweight == sequence > 110% of reference allele length
            if len(seq) / ref_lens[locus] > 1.1:
                status[isolate_id] = "overweight"
                continue

            # Otherwise, it is a complete sequence
            status[isolate_id] = "complete"
        statuses[locus] = pd.Series(status, dtype="category")
    statuses = pd.DataFrame(statuses, dtype="category").sort_index()

    import ipdb; ipdb.set_trace()
    pass

# Toxin gene processing
# Create status df
# Need look-up of all toxin reference alleles
# For each toxin gene file (xxxAxxx.fas - check reference allele is available)
# Open file,
# if no sequence, nd
# if sequence < 80% or > 120% len(allele1), len diff
# if no start codon, or stop codon, or internal stop codons, incomplete CDS
    
if __name__ == '__main__':
    logging.basicConfig(level=logging.DEBUG)
    main()
