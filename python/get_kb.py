#!/usr/bin/env python

import numpy as np
from plot_hist import hist_1d
import matplotlib.pyplot as plt
from scipy.stats import entropy
import sys

if len(sys.argv)<2:
    print "Usage: get_kb.y <DATA1> <DATA2>"
    print "Note: If <DATA2> is not given, <DATA1> will be split in half."
    exit(1)

def calc_kld(p1, p2):
    
    """
    Calculates Kullback-Leibler Divergence.
    """

    valids = np.where((p1>0.001)*(p2>0.001))
    div    = entropy(p1[valids], p2[valids])
    
    return div

def series_kld(data1, data2):

    all_data = np.concatenate((data1, data2))

    _max = np.max(all_data)
    _min = np.min(all_data)

    kde  = hist_1d(all_data,
        xlim=(_min, _max))

    L    = all_data.shape[0]

    series = np.zeros(L-2)
    xx     = np.zeros(L-2)

    for i in range(0,L-2):
        t    = i+2
        kde1 = hist_1d(all_data[:t],
            xlim=(_min, _max))
        series[i] = calc_kld(kde.f, kde1.f)
        xx[i]     = t

    return xx, series

if len(sys.argv)>2:
    file1 = sys.argv[-2]
    file1_data = np.loadtxt(file1)[:,-1]
    file2 = sys.argv[-1]
    file2_data = np.loadtxt(file2)[:,-1]
else:
    file1 = sys.argv[-1]
    file1_data = np.loadtxt(file1)
    L2         = file1_data.shape[0]/2
    file2_data = file1_data[L2:,-1]
    file1_data = file1_data[:L2,-1]

max1 = np.max(file1_data)
max2 = np.max(file2_data)
min1 = np.min(file1_data)
min2 = np.min(file2_data)
        
kde      = hist_1d(np.concatenate((file1_data, file2_data)),
    xlim=(np.min([min1, min2]), np.max([max1, max2])))
kde1 = hist_1d(file1_data,
    xlim=(np.min([min1, min2]), np.max([max1, max2])))
kde2 = hist_1d(file2_data,
    xlim=(np.min([min1, min2]), np.max([max1, max2])))

plt.plot(kde.xx,  kde.f,  label="full dataset")
plt.plot(kde1.xx, kde1.f, label="1st half dataset")
plt.plot(kde2.xx, kde2.f, label="2nd half dataset")
plt.xlabel(r'RMSD $[\AA]$')
plt.ylabel('Density')
plt.legend()
plt.savefig("kde.png", dpi=1000)
plt.clf()

div    = calc_kld(kde1.f, kde2.f)
print "KB divergence D(<DATA1>||<DATA2>): %s" %div

xx, series = series_kld(file1_data, file2_data)
plt.plot(xx, series)
plt.xlabel("No. of Frames")
plt.ylabel("KL-Divergence")
plt.savefig("kld-vs-t.png", dpi=1000)
plt.clf()
np.savetxt("kld-vs-t.dat", np.stack((xx, series), axis=1))