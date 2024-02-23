import sys
import subprocess

class snpdif(object):

    def __init__(self,args):
        self.output = args.output
        self.fasta = str(self.output + 'snps.fasta')

    def calc_snpdif(self):
        parts = []
        f = open(self.output+'script.R','w')
        parts.append('library(Biostrings)')
        parts.append('dna = readDNAStringSet("'+self.fasta+'")')
        parts.append('mat = stringDist(dna, method="hamming")')
        parts.append('mat = as.matrix(mat)')
        parts.append('write.table(mat,file="'+self.output+'snpdif.csv",sep=",")')
        f.write(str('\n'.join(parts)))
        f.close()

        try:
            process = subprocess.Popen(['R', 'CMD', 'BATCH', self.output+'script.R'])
            process.wait()
        except:
            IOError('Could not build R distance matrix')
            sys.exit(1)


