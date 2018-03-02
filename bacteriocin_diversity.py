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
import seaborn as sns
from matplotlib import pyplot as plt
import bidict

def clusters(df):
	clusters = collections.defaultdict(list)
	for column in df.columns:
		if column.startswith("mj_"):
			clusters[column.replace("mj_", "")[:3]].append(column)

	return clusters

def read_export(allele_export):
	df = pd.read_excel(allele_export)
	loci = [locus for loci in clusters(df).values() for locus in loci]
	missing_data = df[df[loci].isnull().any(axis="columns")]["id"]
	if not missing_data.empty:
		logging.warn("Dropped uncurated isolate(s) with id(s): %s", ", ".join(missing_data.astype(str)))
		df.drop(missing_data.index, inplace=True)

	return df


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
@click.option("--output_directory", 
	default="outdir",
	prompt=True,
	type=click.Path(file_okay=False, dir_okay=True, readable=True, writable=True, resolve_path=True))
@click.option('--overwrite', is_flag=True)
def main(allele_export, output_directory, overwrite):
	"""Short description of what the script does

	Long description of what the script does"""
	outdir = prepare_outdir(output_directory, overwrite)

	df = read_export(allele_export)

	c = clusters(df)
	logging.debug("Clusters: %s", c)

	# Summary allele presence across clusters
	summary = pd.DataFrame(index=df.index)
	for cluster, loci in c.items():
		summary[cluster] = (df[loci] != 0).sum(axis="columns") / len(loci)
	summary.to_csv(outdir / "presence_absence_summary.csv")

	# Plot average presence of loci in each cluster
	plt.gcf().set_size_inches(18.5, 10.5)
	sns.heatmap(summary[summary.mean().sort_values().index], yticklabels=False, cmap='YlGnBu')
	plt.savefig(str(outdir / 'presence_absence_heatmap.png'), bbox_inches='tight')
	
	# Replace [s][i] and [s] with EOC and I, respectively
	loci = [locus for loci in c.values() for locus in loci]
	df.loc[:, loci] = df.loc[:, loci].replace(["[S][I]", "[S]"], ["EOC", "I"])
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
	import ipdb; ipdb.set_trace()


if __name__ == '__main__':
	logging.basicConfig(level=logging.DEBUG)
	main()
