import subprocess
import sys
import time

class phylo(object):

    def __init__(self,args):
        self.output = args.output
        self.bootstrap = args.bootstrap

    def phylo_tree(self):
        parts = []
        parts.append('mafft')
        parts.append('--auto')
        parts.append('--quiet')
        parts.append(self.output+'snps.fasta')
        parts.append('>')
        parts.append(self.output+'aligned.fasta')
        mcmd = str(' '.join(parts))

        try:
            print('\n')
            sys.stdout.write('running multiple sequence alignment')
            sys.stdout.flush()
            print('\n')
            process = subprocess.Popen(mcmd,shell=True)
            process.wait()
        except:
            IOError('Could not run mafft: check fasta files')
            sys.exit(1)

        time.sleep(1)

        parts = []
        parts.append('iqtree')
        parts.append('-nt')
        parts.append('AUTO')
        parts.append('-m')
        parts.append('MFP')
        parts.append('-quiet')
        parts.append('-s')
        parts.append(self.output+'aligned.fasta')
        if self.bootstrap == True:
            parts.append('-b')
            parts.append('100')
            icmd = str(' '.join(parts))
            
        try:
            print('\n')
            sys.stdout.write('building phylogenetic tree')
            sys.stdout.flush()
            print('\n')
            process = subprocess.Popen(icmd,shell=True)
            process.wait()
        except:
            IOError('Could not run iqtree: check fasta files')
            sys.exit(1)
