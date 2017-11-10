from __future__ import print_function
import sys
import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import argparse



def plot_rmsd(data,  state):

    fig, ax =  plt.subplots()

    ax.scatter(data['#Frame'], data['RMSD_00000'], alpha=.2)

    data['acc_rmsd'] = data['RMSD_00000'].expanding(2).mean()
    ax.plot(data['#Frame'], data['acc_rmsd'], "r")
    #y_av = movingaverage(dvdl, 200)
    #ax.plot(range(len(dvdl)), y_av,"r")
    ax.set_ylim(0,4)
    ax.set_title(state)
    fig.savefig('rmsd.png')
    plt.close(fig)

def plot_rmsf(data,  state):

    fig, ax =  plt.subplots()
    ax.plot(data.index,data['RMSF'],alpha=.6)
    ax.set_ylim(0,6)
    #y_av = movingaverage(dvdl, 200)
    #ax.plot(range(len(dvdl)), y_av,"r")
    ax.set_title(state)
    fig.savefig('rmsf.png')
    plt.close(fig)


# Plot RMSD
rmsd_file = 'rmsd_bb.1.dat'
rmsf_file = 'rmsf_bb.1.dat'

if os.path.isfile(rmsd_file):
    data = pd.read_csv(rmsd_file, delim_whitespace=True)
    s = os.path.basename(os.path.normpath('.'))
    plot_rmsd(data, s)

if os.path.isfile(rmsf_file):
    data = pd.read_csv(rmsf_file, delim_whitespace=True,names=['AA','RMSF'], skiprows=1)
    s = os.path.basename(os.path.normpath('.'))
    plot_rmsf(data, s)
