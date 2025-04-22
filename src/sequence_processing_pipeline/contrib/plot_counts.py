import matplotlib.pyplot as plt
import re
import sys
import os
import pandas as pd

ex = re.compile(r'_I1_(C5\d\d).fastq.gz.corrected.err_barcode_removed.fastq')

# remove total line from wc
data = [x.strip().split(' ') for x in open(sys.argv[1])][:-1]
plotdata = [(ex.search(i).groups()[0], int(v) / 4) for v, i in data]
sheetdata = dict(plotdata)

ordered = sorted(plotdata, key=lambda x: x[1])
f = plt.figure(figsize=(16, 8))
plt.bar([i for i, _ in ordered], [v for _, v in ordered])
plt.ylabel('I1 reads')
plt.xticks(list(range(len(ordered))), [i for i, _ in ordered], rotation=90)
plt.savefig(sys.argv[3] + '/counts.pdf')

sheet = pd.read_csv(sys.argv[2], dtype=str)
sheet = sheet[~sheet['Lane'].isnull()]
sheet['read_counts'] = [sheetdata[i] for i in sheet['Barcode_ID']]
name = os.path.basename(sys.argv[2]).rsplit('.', 1)[0]
newname = name + '.read_counts.tsv'

sheet.to_csv(sys.argv[3] + '/' + newname, sep='\t', index=False, header=True)
