import csv
from collections import defaultdict

class drugs(object):

	def __init__(self,args,samples,lineage_names,lineage_numbers):
		self.input = args.input
		self.output = args.output
		self.frequency = float(args.frequency)
		self.quality = float(args.quality)
		self.coverage = float(args.coverage)

		self.samples = samples
		self.lineage_names = lineage_names
		self.lineage_numbers = lineage_numbers
		
		self.samples_variants = defaultdict(list)

	def drugger(self):
		print('.... testing drug resistance against database bin/db/drdb.txt ....')
		for sample in self.samples:
			fin = open(self.input + 'snp/' + sample + '.snp','rb')
			reader = csv.DictReader(fin,delimiter = ',')
			for line in reader:
				if float(line['Frequency']) >= self.frequency and float(line['Average quality']) >= self.quality and float(line['Coverage']) >= self.coverage:
					self.samples_variants[sample].append(line['Reference Position'])

		print(self.samples_variants)
