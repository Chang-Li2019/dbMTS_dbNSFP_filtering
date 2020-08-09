
def b2s(byte_obj):
    #byte to string
    return byte_obj.decode('utf-8')



mRNA_dict = {}

# read candidate genes/mRNAs
with open('WWNC.txt','r') as f:
    for i in f:
        i = i.strip('\n')
        mRNA_dict[i] = 0


import gzip, time
Outfile = open('nonsynonymous_SNPs','w')
l = ['carcinoma','sarcoma','leukemia','tumor','cancer','lymphoma','myeloma']
# list of all chromosome names 1-22, M, X, Y
# dbNSFP4.0a_variant.chr1.gz
# dbNSFP4.0a_variant.chr2.gz ... dbNSFP4.0a_variant.chrY.gz
chromes = list(range(1,23))+['M','X','Y']
n = 0
for i in chromes:
    time_now = time.time()
    path = 'dbNSFP4.0b1a_variant.chr{}.gz'.format(i)
    file = gzip.open(path,'r')
    for line in file:
        line = b2s(line)
        line = line.strip('\n').split('\t')
        if line[12] in mRNA_dict:
			if line[112] =='D':
				if line[156]!='.' and line[172]!='.' and float(line[156]) < 0.05 and float(line[172]) < 0.05:
					tag = 0
					for y in l:
						if y in line[364]:# or y in line[366]:
							tag = 1
							break
					if tag == 1:
						Outfile.write('\t'.join(line)+'\n')