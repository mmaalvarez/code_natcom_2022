from random import shuffle


def shuffle_dict(dictionary):
    keys = list(dictionary.keys())
    shuffle(keys)
    shuffled_dict = dict(zip(keys, dictionary.values()))
    return shuffled_dict


ITER = 1000


with open("file_list.txt", "r") as filenames:
	file_list = filenames.read().splitlines()

conditions = ["p53KO_t9_woDOX", \
				"p53KO_t9_DOX", \
				"p53KO_t12_woDOX", \
				"p53KO_t12_DOX", \
				"p53KO_t15_woDOX", \
				"p53KO_t15_DOX", \
				"p53wt_t9_woDOX", \
				"p53wt_t9_DOX", \
				"p53wt_t12_woDOX", \
				"p53wt_t12_DOX", \
				"p53wt_t15_woDOX", \
				"p53wt_t15_DOX"]

files_dict = dict(zip(conditions, file_list))


# original bottomcenter genes per file
genes_dict = dict()

for condition in conditions:
	path = files_dict[condition]

	with open(path, "r") as flute_res:
		for line in flute_res:
			split_line = line.strip().split()
			gene = split_line[0]
			if gene != 'Gene':
				if "bottomcenter" in split_line:
					group = 1
				else:
					group = 0

				try:
					genes_dict[condition][gene] = group
				except KeyError:
					genes_dict[condition] = {gene: group}

genes_list = list(genes_dict["p53KO_t9_woDOX"].keys())


# shuffling
list_n_false_hits_per_iter = []
for n in range(ITER):

	# initialize gene_bottomcenter_counts
	gene_bottomcenter_counts = {"p53KO":dict(), \
								"p53wt":dict()}
	for gene in genes_list:
		gene_bottomcenter_counts["p53KO"][gene] = 0
		gene_bottomcenter_counts["p53wt"][gene] = 0

	shuffled_dict = shuffle_dict(genes_dict)
	for condition in conditions:
		if "KO" in condition:
			for gene in shuffled_dict[condition]:
				gene_bottomcenter_counts["p53KO"][gene] += shuffled_dict[condition][gene]
	for condition in conditions:
		if "wt" in condition:
			for gene in shuffled_dict[condition]:
				gene_bottomcenter_counts["p53wt"][gene] += shuffled_dict[condition][gene]

	# how many genes are 4 in one p53 status and 0 in the other p53 status
	n_false_hits = 0
	for gene in genes_list:
		if gene_bottomcenter_counts["p53KO"][gene] >= 4 and gene_bottomcenter_counts["p53wt"][gene] == 0:
			n_false_hits += 1
		# elif gene_bottomcenter_counts["p53wt"][gene] >= 4 and gene_bottomcenter_counts["p53KO"][gene] == 0:
		# 	n_false_hits += 1
	list_n_false_hits_per_iter.append(n_false_hits)

print(list_n_false_hits_per_iter)

# output
with open("list_n_false_hits_per_iter.tsv", "w") as out:
	out.write(f'{"	".join([str(x) for x in list_n_false_hits_per_iter])}\n')
