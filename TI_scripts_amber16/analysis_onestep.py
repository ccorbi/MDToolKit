from __future__ import print_function
import sys
import os
from collections import defaultdict
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import glob
import re
import argparse

from scipy.integrate import trapz, simps


def plot_rmsd(data, clambda, state):

    fig, ax =  plt.subplots()

    ax.scatter(data['#Frame'], data['RMSD_00000'], alpha=.2)
    #y_av = movingaverage(dvdl, 200)
    #ax.plot(range(len(dvdl)), y_av,"r")
    ax.set_title(clambda+ ' '+ state)
    fig.savefig('./analysis/rmsd/rmsd_{}_{}.png'.format(state,clambda ))


def plot_lambda(dvdl, clambda, state):

    fig, ax =  plt.subplots()

    df = pd.DataFrame(dvdl, columns=['dvdl'])
    df['acc_dvdl'] = df.expanding(2).mean()

    ax.scatter(range(df.shape[0]), df['dvdl'], alpha=.2)

    ax.plot(range(df.shape[0]), df['acc_dvdl'], "r")
    ax.set_title(clambda)
    fig.savefig('./analysis/lambdas/{}_{}.png'.format(state,clambda ))
    plt.close(fig)


def movingaverage(interval, window_size):
    window= np.ones(int(window_size))/float(window_size)
    return np.convolve(interval, window, 'same')


# parse energy output
def parse_energyOut(file_name):
    dVdl = list()
    with open(file_name, 'r') as en_file:
        for line in en_file:
          if  line.startswith('L9') and not 'dV/dlambda' in line:
             dVdl.append(float(line.split()[5]) )

    return dVdl

# Get all lambdas

def islambdafolder(folder):

    if re.match(r'[0-1]\.[0-9][0-9][0-9]',folder):
        return True
    else:
        return False


# print deltas

def get_integration(df):

    df.sort_values('Lambda', inplace=True)
    integration =  trapz(df['DVDL'].get_values(), df['Lambda'].get_values())

    return integration


def parse_TIout(folders,labels ,preprocess_dvdl, SKIP, LIMIT):

    # Parse TI output
    # p_H_data = list()
    # for state in STATES:
    #     for step in STEPS:
    #         folders = glob.glob('./'+state+'/'+step+'/*')
    #labels is a list, wih state and step or only state

    raw_dvdl = defaultdict(list)
    # for lamdba
    for clambda in folders:
        l = os.path.basename(os.path.normpath(clambda))
        if islambdafolder(l):
            energy_files = glob.glob(clambda+'/ti*.en')
            energy_files.sort()

            if len(energy_files):
                for e in energy_files:
                    if e:
                        dvdl = parse_energyOut(e)
                        raw_dvdl[l].extend(dvdl)

                if len(raw_dvdl[l]) >SKIP:
                    plot_lambda(raw_dvdl[l], l, '_'.join(labels))
                    data = np.asarray(raw_dvdl[l])
                    preprocess_dvdl.append(labels+[float(l),data[SKIP:LIMIT].mean(), data[SKIP:LIMIT].std()])


def create_folder(folder):
    '''This only to keep compability with py27, no if_exist'''

    try:
        os.makedirs(folder)
    except:
        pass

################
################
# T O D O
########
# refactor code
# add options RMSD, etc ...
# remove warinings
# return the ddG
# skip lambdas
# add option to skip N time from the calculation of the average dvdl
# add option to setup a maxmim
# add average lamda in the plots
# save data in csv format
# improve RMSD generation, may be use pjtraj to get values
# split plots in folders RMSD, LAmbdas, integration


def get_args():
    ## todo simplify mutations to X#X
    parser = argparse.ArgumentParser()
    parser.add_argument("--extrapol", type=bool, default=True)
    parser.add_argument("--ignore", type=str, default=.1 )
    parser.add_argument("--skip", type=int,default=500 )
    parser.add_argument("--limit", default=None )
    parser.add_argument("--psteps", type=int,default=2500000 )
    #parser.set_defaults(nice=False)
    args = parser.parse_args()
    return args



if __name__ == '__main__':
    args = get_args()

    SKIP = args.skip
    if args.limit:
        LIMIT = int(args.limit)
    else:
        LIMIT = args.limit

    ## CONSTANTS
    STATES = ['complex', 'ligand']
    #STEPS = ['decharge','vdw_bonded', 'recharge']

    # make folder

    create_folder('./analysis')
    create_folder('./analysis/lambdas')
    create_folder('./analysis/rmsd')

    # Parse TI output
    preprocess_dvdl = list()
    for state in STATES:
        #for step in STEPS:

            folders = glob.glob('./'+state+'/*')
            labels=[state]
            parse_TIout(folders, labels , preprocess_dvdl, SKIP, LIMIT)


    df = pd.DataFrame(preprocess_dvdl, columns=['state' ,'Lambda','DVDL', 'rms'])


    # add lambda 0 and 1, lineal extrapolation
    if args.extrapol:
        for state in STATES:


                d = df[(df['state']==state)]

                d.sort_values('Lambda', inplace=True)
                if d[d['Lambda'] == 0.0].empty:

                    l0 = ((d.iloc[0]['Lambda']*d.iloc[1]['DVDL'])-(d.iloc[1]['Lambda']*d.iloc[0]['DVDL']))/(d.iloc[0]['Lambda']-d.iloc[1]['Lambda'])
                    t = pd.DataFrame([[state,0.00,l0,0.00]], columns=['state','Lambda','DVDL', 'rms'])
                    df = df.append(t)

                if d[d['Lambda'] == 1.0].empty:
                    x = ((d.iloc[-2]['Lambda']-1.0)*d.iloc[-1]['DVDL'])+((1.0-d.iloc[-1]['Lambda'])*d.iloc[-2]['DVDL'])
                    y = (d.iloc[-2]['Lambda']-d.iloc[-1]['Lambda'])
                    t = pd.DataFrame([[state,1.00,x/y,0.00]], columns=['state', 'Lambda','DVDL', 'rms'])
                    df = df.append(t)

    # save lambdas
    df.to_csv('./analysis/dvdldata.csv', index=False)

    final = dict()
    fig, ax = plt.subplots(nrows=2, ncols=1, figsize=(6,10))
    #
    #
    delta = open('./analysis/delta.out', 'w')
    for idx,state in enumerate(STATES):
        i = list()
        #for idj, step in enumerate(STEPS):

        report = df[(df['state']==state)]
        integration = get_integration(report)
        #print('# {} \t\t {:.6f}'.format(state, integration))
        i.append(integration)
        report.sort_values('Lambda', inplace=True)
        ax[idx].errorbar(report['Lambda'], report['DVDL'], yerr=report['rms'])
        ax[idx].set_title(state)


        print('dG sum for {}:  \t\t {:.6f}'.format(state,sum(i)))
        print('dG sum for {}:  \t\t {:.6f}'.format(state,sum(i)), file=delta)
        print('-'*30)
        print('\n')
        final[state] = sum(i)

    print('#'*30)
    # ΔGcomplex - ΔGligands
    print('final ddG binding: \t\t {:.6f}'.format( final['complex']-final['ligand']))
    print('final ddG binding: \t\t {:.6f}'.format( final['complex']-final['ligand']), file= delta)

    fig.savefig('./analysis/plot_{}.png'.format('integration'))

    plt.close(fig)
    delta.close()




    # Plot RMSD
    for state in STATES:
        folders = glob.glob('./'+state+'/*')
        raw_H_data = defaultdict(list)

        for clambda in folders:
            if os.path.isfile(clambda+'/rmsd_bb.1.dat'):
                l = pd.read_csv(clambda+'/rmsd_bb.1.dat', delim_whitespace=True)
                s = os.path.basename(os.path.normpath(clambda))
                plot_rmsd(l,s, state)
