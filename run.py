import sys
import os
import argparse

sys.path.append('lib')
from lib import main
from snpdif import snpdif
from lineage import lineage
from phylo import phylo
from metadata import metadata
from itol_parser import itol
from drug_resistance import drugs

def clear(outdir):
	try:
		os.system('rm '+outdir+'*')
	except:
		raise IOError('Could not clear output directory - check permissions')
		sys.exit(1)

def query_genes():
	with open('bin/query_genes.txt') as f:
		metagenes = f.read().splitlines()
	f.close()
	return metagenes

if __name__ == '__main__':

	parser = argparse.ArgumentParser(description = 'FROM CLC VARIANT CALLING TO ML TREE')
	parser.add_argument('-i', '--input', type=str, metavar='', default="input/", help="input directory of SNP files (def = input)")
	parser.add_argument('-o', '--output', type=str, metavar='', default="output/", help="output directory (def = output)")
	parser.add_argument('-f', '--frequency', type=int, metavar='', default=90, help="minimum frequency of SNPs (def = 90)")
	parser.add_argument('-q', '--quality', type=int, metavar='', default=20, help="minimum quality of SNPs (def = 20)")
	parser.add_argument('-c', '--coverage', type=int, metavar='', default=0, help='average coverage of sample to be used (def = 0)')
	parser.add_argument('--ares',dest='ares',default=False, action='store_false',help='exclude antibiotic resistance SNPs (def = False)')
	parser.add_argument('--exon', dest='exon', default=False, action='store_true', help='use only SNPs within exon regions (def = False)')
	parser.add_argument('--metadata', dest='metadata', default=False, action='store_true', help='include metadata file with phylomarkers and specified genes (def = False)')
	parser.add_argument('--phylo', dest='phylo', default=False, action='store_true', help='include metadata file with phylomarkers and specified genes (def = False)')
	parser.add_argument('--bootstrap', dest='bootstrap', default=False, action='store_true', help='allow non-parametric bootstrap of 100 replicates (def = False)')

	args = parser.parse_args()

	clear(args.output)

	imain = main(args)
	imain.snpper()

	ilineage = lineage(args,imain.get_sample_lineages())
	ilineage.lineage()

	isnpdif = snpdif(args)
	isnpdif.calc_snpdif()
	
	if args.metadata == True:
		imetadata = metadata(args,imain.get_samples(),query_genes(),ilineage.get_lineage_names(),ilineage.get_lineage_numbers())
		imetadata.metadata()
		iitol_parser = itol(args,query_genes())
		iitol_parser.annotations()
		#idrugs = drugs(args,imain.get_samples(),ilineage.get_lineage_names(),ilineage.get_lineage_numbers())
		#idrugs.drugger()
	
	if args.phylo == True:
		iphylo = phylo(args)
		iphylo.phylo_tree()

