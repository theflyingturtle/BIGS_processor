#!/usr/bin/env python3

# 			i) allele export (select [s][i])
# 				Counts
# 					Complete - rows where all cells in [loci] =! 0
# 					Partial - rows where >1 cells in [loci] == 0, but not all
# 					Not detected - rows where all cells in [loci] == 0
# 				Presentation
# 					CSV ready for prettifying/R
# 				Final check
# 					Total = 100%

import pandas as pd
import xlrd
import click
import collections
import logging
import os
import shutil
import pathlib
import sys

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
	os.mkdir(output_directory)


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
	prepare_outdir(output_directory, overwrite)

	df = read_export(allele_export)

	c = clusters(df)
	logging.debug("Clusters: %s", c)

	summary = pd.DataFrame(index=df.index)
	for cluster, loci in c.items():
		summary[cluster] = (df[loci] != 0).sum(axis="columns") / len(loci)




	
	import ipdb; ipdb.set_trace()	


if __name__ == '__main__':
	logging.basicConfig(level=logging.DEBUG)
	main()
