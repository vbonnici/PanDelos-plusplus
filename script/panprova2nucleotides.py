import sys

ifasta = sys.argv[1]
igff = sys.argv[2]

ofile = sys.argv[3]

genome = ""
for line in open(ifasta,'r'):
    if line[0] != '>':
        genome += line.strip()


off = open(ofile,'a')

rcmap = {'a':'t', 't':'a', 'c':'g', 'g':'c', 'n':'n', 
    'A':'T', 'T':'A', 'C':'G', 'G':'C', 'N':'N'}
def maprc(c):
    if c in rcmap:
        return rcmap[c]
    else:
        return c

geneid = 0
gene_counts = dict()

for line in open(igff,'r'):
    if line[0] != '#':
        cc = line.strip().split('\t')
        
        #print(cc[3],cc[4],cc[6])

        genome_id = cc[0]
        family = cc[8].split('name=')[1]
        gene_counts[family] = gene_counts.get(family, -1) +1
        gc = gene_counts[family]

        
        start = int(cc[3])
        end = int(cc[4])
        if start > end:
            t = start
            start = end
            end = t
        seq = genome[start:end].upper()
        if cc[6] == '-':
            rc = ""
            for i in range(len(seq)-1, -1,-1):
                rc += maprc(seq[i])
            seq = rc

        off.write(genome_id +'\t' + ( genome_id+':'+str(geneid)+':'+str(gc) ) +'\t'+family+'\n')
        off.write(seq+'\n')

        geneid += 1

off.flush()
off.close()        

