import csv

class exclude(object):

	def __init__(self,args):
		self.ares = args.ares
		self.dr = []
		f = open('bin/db/drdb.txt','r')
		reader = csv.DictReader(f,delimiter='\t')
		for line in reader:
			if '/' in line['Reference_Position']:
				player = line['Reference_Position'].split('/')
				for item in player:
					self.dr.append(item)
			else:				
				self.dr.append(line['Reference_Position'])
		f.close()

		self.repeats = []
		f = open('bin/db/repeats.txt','r')
		reader = csv.DictReader(f,delimiter='\t')
		for line in reader:
			self.repeats.append(line['Start']+':'+line['End'])
		f.close()

		f = open('bin/db/ne_genes.txt','r')
		reader = csv.DictReader(f,delimiter='\t')
		for line in reader:
			self.repeats.append(line['Start']+':'+line['End'])
		f.close()		
	
	def get_dr(self):
		if self.ares == True:
			return self.dr
		else:
			return ''

	def test_repeat(self,snp):
		self.snp = int(snp)
		for k in self.repeats:
			s = int(k.split(':')[0])
			e = int(k.split(':')[1])
			if self.snp >= s and self.snp <= e:
				player = True
				break
			else:
				player = False
		return player

	def proximity(self,snps):
		d = []
		for x in snps:
			for y in snps:
				if 0 < int(int(y) -int(x)) <= 50:
					d.append(x)
					d.append(y)
		d = set(d)
		return d
			
			
		
