#!/usr/bin/env python3

import matplotlib
matplotlib.use("Agg")

import pandas as pd
import xlrd
import click
import collections
import logging
import os
import shutil
import pathlib
import sys
import bidict
import seaborn as sns

from collections import defaultdict
from matplotlib import pyplot as plt

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio.Data.CodonTable import TranslationError


def clusters(df):
    clusters = collections.defaultdict(list)
    for column in df.columns:
        if column.startswith("mj_"):
            clusters[column.replace("mj_", "")[:3]].append(column)

    return clusters


def read_export(allele_export):
    df = pd.read_excel(allele_export, dtype=str)
    loci = [locus for loci in clusters(df).values() for locus in loci]

    # Replace "[S][I]" and "[S]" with "EOC" (end of contig) and "Q" (questionable), respectively
    df.loc[:, loci] = df.loc[:, loci].replace(["[S][I]", "[S]"], ["EOC", "Q"])

    # Drop rows where loci appear to be uncurated
    # Must check =="nan" since we read everything as strings, which we did for the string replacement above
    missing_data = df[(df[loci] == "nan").any(axis="columns")]["id"]
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


def read_toxins(toxin_seqs):
    toxins = {}
    for f in pathlib.Path(toxin_seqs).glob("*.fas"):
        dedup_nucs = defaultdict(list)
        for record in SeqIO.parse(str(f), "fasta"):
            dedup_nucs[str(record.seq)].append(record.id)
        toxins[f.stem] = dedup_nucs

    return toxins


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
# TODO: don't import toxins again, instead identify them in allele_seqdir
@click.argument("toxin_seqdir", type=click.Path(exists=True, file_okay=False))
@click.option("--output_directory", 
    default="outdir",
    prompt=True,
    type=click.Path(file_okay=False, dir_okay=True, readable=True, writable=True, resolve_path=True))
@click.option('--overwrite', is_flag=True)
def main(allele_export, ref_seqs, allele_seqdir, toxin_seqdir, output_directory, overwrite):
    """TO DO... Short description of what the script does

    Long description of what the script does"""
    outdir = prepare_outdir(output_directory, overwrite)

    df = read_export(allele_export)

    c = clusters(df)
    logging.debug("Clusters: %s", c)

    # Summary allele presence across clusters
    cluster_summary = pd.DataFrame(index=df.index)    
    for cluster, loci in c.items():
        cluster_summary[cluster] = (df[loci] != "0").sum(axis="columns") / len(loci)
    cluster_summary.to_csv(outdir / "cluster_presence_absence_summary.csv")

    # Plot average presence of loci in each cluster
    plt.gcf().set_size_inches(18.5, 10.5)
    sns.heatmap(cluster_summary[cluster_summary.mean().sort_values().index], yticklabels=False, cmap='YlGnBu')
    plt.savefig(str(outdir / 'presence_absence_heatmap.png'), bbox_inches='tight')

    # Summary allele presence within clusters
    gene_list = []
    for cluster, loci in c.items():
          gene = (df[loci] != "0").sum(axis="rows") / len(df.index)
          gene_list.append(gene)
    gene_list = pd.concat(gene_list)
    gene_prev = pd.DataFrame({'locus':gene_list.index, 'proportion_present':gene_list.values})
    gene_prev['cluster'] = gene_prev['locus'].replace('mj_', '', regex=True).str[:3]
    gene_prev = gene_prev.reindex(sorted(gene_prev.columns), axis=1)
    gene_prev.to_csv(outdir / "gene_prevalence_summary.csv")

    # Summary of bacteriocin profiles (as for MLST)
    bact_STs = {}
    isolate_STs = {}
    for cluster, loci in c.items():
        # This whole hacky section is because we need to assign the profile '0-0-...-0' the key -9
        profiles = df.set_index("id").loc[:, loci].astype(str).apply("-".join, axis=1)
        is_zero = lambda profile: all(allele == '0' for allele in profile.split("-"))
        has_an_all_zero_profile = any(is_zero(profile) for profile in profiles.unique())
        nonzero_profiles = (profile for profile in profiles.unique() if not is_zero(profile))

        bact_STs[cluster] = bidict.bidict(enumerate(nonzero_profiles, 1))
        if has_an_all_zero_profile:
            bact_STs[cluster][-9] = next(profile for profile in profiles.unique() if is_zero(profile))
        isolate_STs[cluster] = profiles.map(bact_STs[cluster].inv.get).astype(int)

    pd.DataFrame(isolate_STs).to_csv(outdir / "isolate_STs.csv")
    pd.DataFrame.from_dict(bact_STs, orient='index').to_csv(outdir / "bact_STs.csv")


    # Process nucleotide sequences
    ref_lens = read_ref_seqs(ref_seqs)
    isolate_alleles = read_isolate_alleles(allele_seqdir)
    statuses = {}
    for locus, alleles in isolate_alleles.items():
        status = {}
        allele_ids = df.set_index("id")[locus].to_dict()
        for isolate_id, seq in alleles.items():
            # Absent == sequence starts with a gap (-)
            if seq.startswith("-"):
                status[isolate_id] = "absent"
                continue

            # EOC == sequence 
            try:
                if allele_ids[str(isolate_id)] == "EOC":
                    status[isolate_id] = "eoc"
                    continue
            except KeyError:
                logging.warn("Failed to look up isolate %s in allele export", isolate_id)

            # Fragment == sequence < 90% of reference allele length
            if len(seq) / ref_lens[locus] < 0.9:
                status[isolate_id] = "fragment"
                continue

            # Overweight == sequence > 110% of reference allele length
            if len(seq) / ref_lens[locus] > 1.1:
                status[isolate_id] = "atypical_length"
                continue

            # Incomplete == sequence does not translate
            # This deliberately comes after fragment and overweight
            try:
                seq.translate(cds=True)
            except TranslationError:
                status[isolate_id] = "incomplete"
                continue

            # Otherwise, it is a complete sequence
            status[isolate_id] = "complete"
        statuses[locus] = pd.Series(status, dtype="category")
    statuses = pd.DataFrame(statuses, dtype="category").sort_index()

    summary = statuses.apply(pd.Series.value_counts).fillna(0).astype(int)
    summary.to_csv(outdir / "locus_status_summary.csv")

    pass

    # Toxin gene processing
    # Nucleotides
    toxins = read_toxins(toxin_seqdir)
    
    tot_nuc_counts = {}
    for toxin, alleles in toxins.items():
        ind_nuc_counts = {}
        for seq, ids in alleles.items():
            tmp_count = len(ids)
            ind_nuc_counts[seq] = tmp_count
        tot_nuc_counts[toxin] = ind_nuc_counts

    unique_toxin_nucs = {}
    for locus, alleles in toxins.items():
        unique_nucs = {}
        for seq, ids in alleles.items():
            representative_id = ids[0]
            count = tot_nuc_counts[locus][seq]
            unique_nucs[seq] = f"{representative_id}_{count}"
        unique_toxin_nucs[locus] = unique_nucs
        
    # Save unique nucleotide sequences
    # FASTA headers follow format of ["representative id"_"number of isolates with this allele"]
    nuc_dir = outdir / "unique_toxin_nucs"
    nuc_dir.mkdir()
    for locus, unique_nucs in unique_toxin_nucs.items():
        nucs_out = (SeqRecord(Seq(k, IUPAC.IUPACAmbiguousDNA), id=v, description="") for k, v in unique_nucs.items())
        SeqIO.write(nucs_out, str(nuc_dir / f"{locus}_unique_nucs.fas"), "fasta")

    # Amino acids
    unique_toxin_aminos = {}
    for f in pathlib.Path(nuc_dir).glob("*.fas"):
        unique_aminos = defaultdict(list)
        for record in SeqIO.parse(str(f), "fasta"):
            record.seq = record.seq.translate()
            unique_aminos[str(record.seq)].append(record.id)
        unique_toxin_aminos[f.stem] = unique_aminos

    amino_dir = outdir / "unique_toxin_aminos"
    amino_dir.mkdir()
    for locus, unique_aminos in unique_toxin_aminos.items():
        aminos_out = (SeqRecord(Seq(k, IUPAC.ExtendedIUPACProtein), id="|".join(v), description='') for k, v in unique_aminos.items())
        SeqIO.write(aminos_out, str(amino_dir / f"{locus}_unique_aminos.fas"), "fasta")

    import ipdb; ipdb.set_trace()
        
    #     dedup_aminos = defaultdict(list)
    #     for record in SeqIO.parse("nucs.fas", "fasta"):
    #         record.seq = record.seq.translate()
    #         dedup_aminos[str(record.seq)].append(record.id)

    #     aminos_out = (SeqRecord(Seq(k, IUPAC.ExtendedIUPACProtein), id="|".join(v), description='') for k, v in dedup_aminos.items())
    #     SeqIO.write(aminos_out, "aminos.fas", "fasta")
 

if __name__ == '__main__':
    logging.basicConfig(level=logging.DEBUG)
    main()
