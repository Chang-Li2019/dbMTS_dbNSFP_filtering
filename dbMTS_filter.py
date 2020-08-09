from datetime import datetime
import itertools
import gzip
startTime = datetime.now()
# Input files
dbmts = open('dbMTS','r')
mir_file = open('miRNA_dict.txt','r')
cand = open('Bai_mRNA_miRNA.txt','r')
gene_to_trans = open('mart_export.txt','r')

# Output file
Outfile = open('Pairs_annotated','w')


def match_key(str1,str2):
	mir_list = str1.split('|')
	trans_list = str2.split('|')
	print_line = 0
	assert(len(mir_list)==len(trans_list))
	pair_observed = []
	for i in range(len(mir_list)):
		if mir_list[i][-3:] in ['-5p','-3p']:
			mir = mir_list[i][0:-3].lower()
		else:
			continue
		trans = trans_list[i]
		if trans in gene:
			key = gene[trans]+'_'+mir_new
			if key in cand_pairs:
				print_line = 1
				key = key.split('_')
				key = ':'.join(key)
				pair_observed.append(key)
	return print_line, pair_observed

gene = {}
for line in gene_to_trans:
	line = line.strip('\n').split('\t')
	gene[line[0]]= line[1].upper()

cand_pairs = {}
for line in cand:
	line = line.strip('\n').split('\t')
	key = line[0].upper()+'_'+line[1].lower()
	cand_pairs[key] = 0

mir = {}
for line in mir_file:
	line = line.strip('\n').split('\t')
	mir[line[2]] = line[1]

for line in dbmts:
	line = line.strip('\n').split('\t')
	key = line[4]+'_'+line[5]+'_'+line[6]+'_'+line[7]
	if key in hgmd:
		count = 0
		ref_pairs = []
		alt_pairs = []
		if count ==0:
			count = match_key(line[162],line[163])
		if count ==0:
			count = match_key(line[175],line[176])
		if count ==1:
			if ref_pairs == []:
				alt_pairs = set(alt_pairs)
				alt_pairs= '|'.join(alt_pairs)
				Outfile.write('\t'.join(line)+'\t'+'.'+'\t'+alt_pairs+'\t'+'\t'.join(hgmd[key])+'\n')
			if alt_pairs == []:
				ref_pairs = set(ref_pairs)
				ref_pairs= '|'.join(ref_pairs)
				Outfile.write('\t'.join(line)+'\t'+ref_pairs+'\t'+'.'+'\t'+'\t'.join(hgmd[key])+'\n')
			if ref_pairs!=[] and alt_pairs != []:
				alt_pairs = set(alt_pairs)
				ref_pairs = set(ref_pairs)
				ref_pairs= '|'.join(ref_pairs)
				alt_pairs= '|'.join(alt_pairs)
				Outfile.write('\t'.join(line)+'\t'+ref_pairs+'\t'+alt_pairs+'\t'+'\t'.join(hgmd[key])+'\n')
