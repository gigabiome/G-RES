import os
import csv
from collections import defaultdict

class itol(object):

	def __init__(self,args,genes):
		self.output = args.output
		self.genes = genes
		self.samples = []
		fmetadata = open(self.output + 'metadata.tab','rb')
		reader = csv.DictReader(fmetadata,delimiter='\t')
		for line in reader:
			self.samples.append(line['Sample'])
		self.samples = set(self.samples)
		
	def annotations(self):
		print('creating itol annotations for ....\n\n')
		for gene in self.genes:
			print(gene + '\n')
			fout = open(self.output + str(gene.replace(' ','_')) + '_itol.txt','w')
			fout.write('DATASET_BINARY\nSEPARATOR TAB\nDATASET_LABEL\t'+str(gene.replace(' ','_'))+'\n')
			fout.write('COLOR\t#ff0000\n\n')
			fout.write('SHOW_LABELS	1\nFIELD_SHAPES\t')
			mutations = []
			sample_mutations = defaultdict(list)
			for sample in self.samples:
				fmetadata = open(self.output + 'metadata.tab','rb')
				reader = csv.DictReader(fmetadata,delimiter='\t')
				for line in reader:
					if str(line['Sample']) == sample and str(line['Name']) == gene:
						aa_change = str(line['Amino acid change']).split('/')
						if len(aa_change) <= 5:
							for x in aa_change:
								mutations.append(x)
								sample_mutations[sample].append(x)
				fmetadata.close()
			mutations = sorted(set(mutations))
			mutations = filter(None,mutations)
			if 'None' in mutations:
				mutations.remove('None')
			mut_len = len(mutations)
			x = 0
			while x < mut_len:
				fout.write('2\t')
				x += 1
			fout.write('\nFIELD_COLORS\t')
			x = 0
			while x < mut_len:
				fout.write('black\t')
				x += 1
			fout.write('\nFIELD_LABELS\t')
			for item in mutations:
				fout.write(item + '\t')
			fout.write('\n\nDATA\n')
			for sample in sample_mutations:
				fout.write(sample + '\t')
				for mutation in mutations:
					if mutation in sample_mutations[sample]:
						fout.write('1\t')
					else:
						fout.write('0\t')
				fout.write('\n')
		print('itol files created .... \n\n')
