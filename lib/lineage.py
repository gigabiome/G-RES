import pandas as pd

class lineage(object):

    def __init__(self,args,sample_lineages):	
        self.input = args.input
        self.output = args.output
        self.sample_lineages = sample_lineages
        self.log = open(self.output + 'log.txt','a')

        df = pd.read_csv('bin/db/lineages.bed', sep='\t', na_values=['?'])	
        self.db = df.set_index('Start', drop=False)

    def lineage(self):
        f = open(self.output + 'lineages.tab', 'w')
        f.write('Sample\tLineage\tLineage_name\tComplete_lineage\n')
        self.log.write('\nThe following sample did not have lineage markers:\n')
        self.lin_names = {}
        self.lin_numbers = {}
        for sample in self.sample_lineages:
            l = {}
            sample_df = pd.read_csv(self.input + 'snp/' + sample + '.snp',sep=',')
            sample_df = sample_df.set_index('Reference Position',drop=False)
            try:
                for lpos in self.sample_lineages[sample]:
                    lpos = int(lpos)
                    if sample_df.loc[lpos,'Reference'] == self.db.loc[lpos,'Ref']:
                        if sample_df.loc[lpos,'Allele'] == self.db.loc[lpos,'Allele']:
                            l[self.db.loc[lpos,'Lineage_number']] = self.db.loc[lpos,'Lineage_name']
            except:
                ValueError
                self.log.write('More than one SNP at a lineage marker: '+str(sample)+'\n')
            dom = ''
            gl = 0
            for k in l:
                x_len = len(k)
                if x_len > gl:
                    gl = x_len
                    dom = k
            try:
                f.write(sample+'\t'+dom+'\t'+str(l[dom])+'\t'+str(l)+'\n')
                self.lin_names[sample] = str(l[dom])
                self.lin_numbers[sample] = str(dom)
            except:
                KeyError
                self.log.write(sample + '\n')
                self.lin_names[sample] = 'N/A'
                self.lin_numbers[sample] = 'N/A'
            f.close()
	
    def get_lineage_names(self):
        return self.lin_names

    def get_lineage_numbers(self):
        return self.lin_numbers
