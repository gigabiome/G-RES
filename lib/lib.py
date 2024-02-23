import os 
import sys
import csv
import subprocess
from collections import defaultdict, OrderedDict
from exclude import exclude

class main(object):

	def __init__(self,args):
		self.input = args.input
		self.output = args.output
		self.frequency = float(args.frequency)
		self.quality = float(args.quality)
		self.coverage = float(args.coverage)
		self.exon = args.exon
		self.ares = args.ares
		self.metadata = args.metadata
		self.log = open(self.output + 'log.txt','a')

		self.indirs = ['snp/','cov/','stat/']

		self.iexclude = exclude(args)
		self.dr = self.iexclude.get_dr()

		self.lineage_dict = defaultdict(list)
		self.lineage_snps = []
		f = open('bin/db/lineages.bed','rb')
		reader = csv.DictReader(f,delimiter='\t')
		for line in reader:
			self.lineage_snps.append(line['Start'])
		f.close()

		#remove special characters
		scmds = []
		for indir in self.indirs:
			scmds.append("rename 's/ /_/g' "+self.input+indir+"*")
			scmds.append("rename 's/\(/_/g' "+self.input+indir+"*")
			scmds.append("rename 's/\)/_/g' "+self.input+indir+"*")

		for scmd in scmds:
			try:
				subprocess.Popen(scmd, shell=True).wait()
			except:
				raise IOError('Could not rename files: '+scmd+'!!!')
				sys.exit(1)

		#rename files
		for indir in self.indirs:
			for file in os.listdir(self.input + indir):
				filehandle = open(self.input+indir+file,'rb')
				try:
					dialect = csv.Sniffer().sniff(filehandle.read(1024))
					parts = file.split('.')
					parts = parts[0].split('_')
					new = str('_'.join(parts[0:2])+'.'+indir.replace('/','')).strip()
					cmd = 'mv '+self.input+indir+file+' '+self.input+indir+new
					try:
						subprocess.Popen(cmd, shell=True).wait()
					except:	
						raise IOError('File '+file+' could not be renamed !!!')
						sys.exit(1)
				except:
					raise IOError('File '+file+' is not a csv file')
					sys.exit(1)

		sfiles = os.listdir(self.input + 'snp/')
		self.snp_files = [self.input +  'snp/' + snp_file for snp_file in sfiles]
		cfiles = os.listdir(self.input + 'cov/')
		self.cov_files = [self.input +  'cov/' + cov_file for cov_file in cfiles]

	def snpper(self):

		self.reference_dict = self.reference()
		fphy = open(self.output + 'snps.phy','w')
		ffasta = open(self.output + 'snps.fasta','w')
		fphy.write(str(len(self.passed_cov_check))+'\t'+str(len(self.reference_dict))+'\n')

		#create alignment file
		print('\n building alignment files ... \n')

		ffasta.write('>NC_000962_H37Rv\n')
		for ref in self.reference_dict:
			ffasta.write(self.reference_dict[ref])
		ffasta.write('\n')
		for snp_file in self.passed_cov_check:
			print(snp_file)
			sample = str(snp_file.split('/')[-1].replace('.snp',''))
			fphy.write(sample + '\t')
			ffasta.write('>' + sample +'\n')

			self.samples_dict = {}
			f = open(snp_file,'rb')
			reader = csv.DictReader(f,delimiter=',')
			for line in reader:
				pos = line['Reference Position']
				vtype = line['Type']
				ref = line['Reference']
				allele = line['Allele']
				cov = float(line['Coverage'])
				freq = float(line['Frequency'])
				qual = float(line['Average quality'])
				anno = line['Overlapping annotations']
				if vtype == 'SNV' and freq >= self.frequency and qual >= self.quality and cov >= self.coverage:
					if self.exon == True:
						if anno:
							self.samples_dict[pos] = allele
					else:
						self.samples_dict[pos] = allele
			f.close()
			for ref_position in self.reference_dict:
				if ref_position in self.samples_dict:
					fphy.write(self.samples_dict[ref_position])
					ffasta.write(self.samples_dict[ref_position])
					if ref_position in self.lineage_snps:
						self.lineage_dict[sample].append(ref_position)
				else:
					fphy.write(self.reference_dict[ref_position])
					ffasta.write(self.reference_dict[ref_position])
			fphy.write('\n')
			ffasta.write('\n')
		fphy.close()
		ffasta.close()

	def reference(self):

		self.passed_cov_check = self.calc_cov()
		self.reference_dict = {}
		self.reference_sorted = OrderedDict()
		print('\n determining global distribution of SNPs ... \n')
		sys.stdout.write('<')
		sys.stdout.flush()
		for sample in self.passed_cov_check:
			sys.stdout.write('-')
			sys.stdout.flush()
			f = open(sample,'rb')
			reader = csv.DictReader(f,delimiter=',')
			for line in reader:
				pos = line['Reference Position']
				vtype = line['Type']
				ref = line['Reference']
				allele = line['Allele']
				cov = float(line['Coverage'])
				freq = float(line['Frequency'])
				qual = float(line['Average quality'])
				anno = line['Overlapping annotations']
				aa = line['Amino acid change']
				rep = self.iexclude.test_repeat(pos)
				if pos not in self.reference_dict:
					if vtype == 'SNV' and freq >= self.frequency and cov >= self.coverage and qual >= self.quality and pos not in self.dr and rep == False:
						if self.exon == True:
							if anno and not aa:
								self.reference_dict[pos] = ref
						else:
							self.reference_dict[pos] = ref
			f.close()
		sys.stdout.write('>')
		sys.stdout.flush()
		prox = list(self.iexclude.proximity(self.reference_dict))
		d_ref = {}
		for key in self.reference_dict:
			if key not in prox:
				d_ref[key] = self.reference_dict[key]
		for key in sorted(d_ref.iterkeys(),key=int):
			self.reference_sorted[key] = d_ref[key]
		return d_ref

	def calc_cov(self):
		self.passed_cov_check = []
		self.log.write('The following samples were not processed:\n')
		for infile in self.snp_files:
			f = open(infile,'rb')
			reader = csv.DictReader(f,delimiter=',')
			self.sample_cov = l = 0
			for line in reader:
				l += 1
				self.sample_cov = self.sample_cov + float(line['Coverage'])
			try:
				self.sample_avg_cov = float(self.sample_cov / l)
				if self.sample_avg_cov >= 10:
					self.passed_cov_check.append(infile)
				else:
					self.log.write(infile + '\n')
			except:
				ZeroDivisionError
				self.log.write(infile + '\n')
		return self.passed_cov_check

	def get_sample_lineages(self):
		return self.lineage_dict

	def get_samples(self):
		samples = []
		for sample in self.passed_cov_check:
			sample = str(sample.split('/')[-1].replace('.snp','')).strip()
			samples.append(sample)
		return samples
