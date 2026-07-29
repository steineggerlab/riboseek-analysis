#! /usr/bin/env python

'''
this script adjusts the bin sizes and takes as input the file that parsed
the results (use get_TP_FP_SCOP.py for scop benchmarking data,
get_FP_shuffle.py for shuffled data and get_evalues_3dblast.py for 3D-BLAST data)

usage:
python3 prep_plot_evalues.py RESULT_FILE DBSIZE TOOL_NAME OUTPUT_FILE

'''

import numpy as np
import sys
import os
import matplotlib.pyplot as plt

result_file = sys.argv[1]
dbsize = int(sys.argv[2])
tool_name  = sys.argv[3]
output_file = sys.argv[4]

FP_evalue = []
scatter_data = []

with open(result_file) as input:
    for line in input:
        data = line.split(',')
        if(data[3].split('\n')[0] == '0'):
            eval = np.float128(data[2])
            # eval = eval ** 0.32
            # FP_evalue.append(eval / 1000)
            FP_evalue.append(eval)

bins = []
maxiFS = max(FP_evalue)
miniFS = min(FP_evalue)
bins.append(0)
stepSize = 1E-50
while stepSize < maxiFS:
    bins.append(stepSize)
    stepSize = stepSize*1.5

hist, bin_edges = np.histogram(FP_evalue,bins)

# Cumulative sum of hist
hist = np.cumsum(hist)

x = []
diag = []

# adjusted bin lengths
for i in range(1,len(bins)):
    x.append(1*(bins[i-1] + (bins[i] - bins[i-1])/2))

fig, ax = plt.subplots()

scatter_y = [item/dbsize for item in hist]

ax.scatter(x, scatter_y, alpha=0.5, color='blue', label=tool_name)
# x and y to 0 to 10
ax.set_xlim([miniFS, 1000])
ax.set_ylim([miniFS, 1000])
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel('Reported e-value')
ax.set_ylabel('Observed false positives per query')
ax.set_title('Reported vs Observed false positives per query (' + tool_name + ')')
# y=x
ax.plot(x, x, color='red', label='y=x')

# Save the figure
plt.savefig(output_file + ".png")

for j in range(0, len(x)):
    scatter_data.append(str(x[j]) + "," + str(scatter_y[j]) + "," + tool_name + "\n")

file = open(output_file, 'w')
for line in scatter_data:
    file.write(line)
file.close()
