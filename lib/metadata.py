import csv,os,sys
import pandas as pd
from collections import defaultdict

class metadata(object):

    def __init__(self,args,samples,genes,lineage_names,lineage_numbers):
        self.input = args.input
        self.output = args.output
        self.frequency = float(args.frequency)
        self.quality = float(args.quality)
        self.coverage = float(args.coverage)

        self.genes = genes
        self.lineage_names = lineage_names
        self.lineage_numbers = lineage_numbers
        with open('bin/coverage_columns.txt') as f:
            self.covcols = sorted(f.read().splitlines())
        f.close()
        with open('bin/variant_columns.txt') as f:
            self.snpcols = sorted(f.read().splitlines())
        f.close()
        with open('bin/statistics_columns.txt') as f:
            self.statcols = sorted(f.read().splitlines())
            f.close()

        self.fmetadata = open(self.output + 'metadata.tab','w')
        self.fmetadata.write('Sample\tLineage_name\tLineage_number\t')
        for k in self.covcols:
            self.fmetadata.write(k+'\t')
        for k in self.statcols:
            self.fmetadata.write(k+'\t')
        for k in self.snpcols:
            self.fmetadata.write(k+'\t')
        self.fmetadata.write('\n')
        
        self.samples = []
        flog = open(args.output+'log.txt','a')
        flog.write('For the following files one of the files (snp / cov / stat) is not available and will not be in metadata file\n')
        for sample in samples:
            fcov = str(self.input+'cov/'+sample+'.cov').strip()
            fsnp = str(self.input+'snp/'+sample+'.snp').strip()
            fstat = str(self.input+'stat/'+sample+'.stat').strip() 
            if os.path.exists(fcov) and os.path.exists(fsnp) and os.path.exists(fstat):
                self.samples.append(sample)
            else:
                flog.write(sample+'\n')
        flog.close()

    def metadata(self):
        print('\n')
        sys.stdout.write('building metadata file for :')
        sys.stdout.flush()
        print('\n')
        for gene in self.genes:
            sys.stdout.write('<----'+gene+'---->')
            sys.stdout.flush()
            print('\n')
            for sample in self.samples:

                snp_dict = defaultdict(list)
                df = pd.read_csv(self.input+'cov/'+sample+'.cov',sep=',')
                df_cov = df.set_index('Name',drop=False)

                self.fmetadata.write(sample +'\t')
                try:
                    self.fmetadata.write(self.lineage_names[sample]+'\t')
                    self.fmetadata.write(self.lineage_numbers[sample]+'\t')
                except KeyError:
                    self.fmetadata.write('N/A\tN/A\t')

                for c in self.covcols:
                    self.fmetadata.write(str(df_cov.loc[gene,c]))
                    self.fmetadata.write('\t')

                f = open(self.input+'stat/'+sample+'.stat','rb')
                reader = csv.DictReader(f,delimiter=',')
                for line in reader:
                    for c in self.statcols:
                        self.fmetadata.write(str(line[c].strip()))
                        self.fmetadata.write('\t')
                f.close()
                f = open(self.input+'snp/'+sample+'.snp','rb')
                reader = csv.DictReader(f,delimiter=',')
                for line in reader:
                    if gene in line['Overlapping annotations']:
                        if float(line['Frequency']) >= self.frequency and float(line['Average quality']) >= self.quality and float(line['Coverage']) >= self.coverage:
                            for c in self.snpcols:
                                if c == 'Average quality':
                                    player = round(float(line[c]),2)
                                    snp_dict[c].append(str(player))
                                elif c == 'Frequency':
                                    player = round(float(line[c]),2)
                                    snp_dict[c].append(str(player))
                                elif c == 'Amino acid change':
                                    player = str(line[c].split('.')[-1]).strip()
                                    snp_dict[c].append(player)
                                elif c == 'Coding region change':
                                    player = str(line[c].split('.')[-1]).strip()
                                    snp_dict[c].append(player)
                                else:
                                    snp_dict[c].append(str(line[c]))
                for c in sorted(snp_dict):
                    self.fmetadata.write(str('/'.join(snp_dict[c])).strip())
                    self.fmetadata.write('\t')
                    self.fmetadata.write('\n')
