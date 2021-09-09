# query = ["dfg","456","AAAAGG"]
barcode_umi = tuple(query[0:2])
sequence = query[2]
if barcode_umi not in list(all.keys()):
	all[barcode_umi] = {sequence : 1}
else:
	if sequence not in list(all[barcode_umi].keys()):
		all[barcode_umi][sequence] = 1
	else:
		all[barcode_umi][sequence] = all[barcode_umi][sequence] + 1

winners = set()
for barcode_umi in all.keys():
	cell = all[barcode_umi]
	seq = [seq for seq,count in cell.items() if count == max(cell.values())] # The key of the max value
	query_winners = tuple([barcode_umi[0], barcode_umi[1], i])
	if len(seq) == 1:
		winners.add(query_winners)
	else:
		for i in seq:
			winners.add(query_winners)
