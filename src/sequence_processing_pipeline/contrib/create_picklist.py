import os
# from metapool.metapool import *
from sys import argv
import pandas as pd
import matplotlib.pyplot as plt
from metapool.metapool import (read_survival, make_2D_array,
                               calculate_iseqnorm_pooling_volumes,
                               format_pooling_echo_pick_list)
import seaborn as sns

input_sheet_filename = argv[1]

plate_df_w_reads = pd.read_csv(input_sheet_filename, sep='\t')
plate_df_w_reads['Blank'] = [True if 'blank' in s.lower() else False
                             for s in plate_df_w_reads['Sample_Name']]
reads_column = 'read_counts'

well_col = 'Sample_Well'
assert reads_column in plate_df_w_reads.columns

f, ((ax1, ax2), (ax3, ax4)) = plt.subplots(nrows=2, ncols=2, figsize=(8, 8))
# evenness plot
rmax = int(round(plate_df_w_reads[reads_column].max(), -2))

foo = read_survival(plate_df_w_reads.loc[plate_df_w_reads['Blank'] is True,
                                         reads_column],
                    label='Blanks',
                    rmax=rmax)

bar = read_survival(plate_df_w_reads.loc[plate_df_w_reads['Blank'] is False,
                                         reads_column],
                    label='Samples',
                    rmax=rmax)

survival_df = pd.concat([foo, bar])

ax3.set_xlabel(reads_column)
ax3.set_ylabel('Samples')
survival_df.plot(color=['coral', 'steelblue'], ax=ax1)
ax1.set_xlabel(reads_column)
ax1.set_ylabel('Samples')

# Histogram
sns.histplot(plate_df_w_reads[reads_column], ax=ax3)

# Boxplot
sns.boxplot(x="Blank", y=reads_column, data=plate_df_w_reads, ax=ax4)
sns.stripplot(x="Blank", y=reads_column, data=plate_df_w_reads, ax=ax4,
              size=3, color='black', alpha=0.5)

plt.tight_layout()
plt.savefig(input_sheet_filename + '.comboplot.pdf')

pdfn = calculate_iseqnorm_pooling_volumes(plate_df_w_reads,
                                          dynamic_range=20,
                                          normalization_column=reads_column)
plt.savefig(input_sheet_filename + '.normalizedplot.pdf')

vols = make_2D_array(pdfn,
                     data_col='iSeq normpool volume',
                     well_col=well_col).astype(float)

# Write the picklist as .csv
picklist_fp = input_sheet_filename + '.picklist.csv'

if os.path.isfile(picklist_fp):
    print("Warning! This file exists already.")

picklist = format_pooling_echo_pick_list(vols, max_vol_per_well=30000)

with open(picklist_fp, 'w') as f:
    f.write(picklist)
