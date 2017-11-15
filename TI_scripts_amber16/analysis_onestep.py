from __future__ import print_function

import argparse
import glob
import os
import re
import sys
import warnings
from collections import defaultdict

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.integrate import simps, trapz

warnings.filterwarnings("ignore")

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


def parse_DVDL_MD(folders,labels ,preprocess_dvdl, skip, limit, ignore, no_plot):

    # Parse TI output
    raw_dvdl = defaultdict(list)
    # for lamdba
    for clambda in folders:
        l = os.path.basename(os.path.normpath(clambda))
        if islambdafolder(l) and l != ignore :
            energy_files = glob.glob(clambda+'/ti*.en')
            energy_files.sort()

            if len(energy_files):
                for e in energy_files:
                    if e:
                        dvdl = parse_energyOut(e)
                        raw_dvdl[l].extend(dvdl)

                if len(raw_dvdl[l]) >skip:
                    if not no_plot:
                        plot_lambda(raw_dvdl[l], l, '_'.join(labels))
                    data = np.asarray(raw_dvdl[l])
                    preprocess_dvdl.append(labels+[float(l),data[skip:limit].mean(), data[skip:limit].std()])


def extrapolation(df, STATES):
    '''If lambda 0 or lambda 1 is missing, add value by extrapolation'''
    for state in STATES:

        d = df[(df['state']==state)]

        d.sort_values('Lambda', inplace=True)
        if d[d['Lambda'] == 0.0].empty:
            
            print("WARNING: Lambda 0 on {}, empty, data point will be interpolated".format(state))
            try:
                l0 = ((d.iloc[0]['Lambda']*d.iloc[1]['DVDL'])-(d.iloc[1]['Lambda']*d.iloc[0]['DVDL']))/(d.iloc[0]['Lambda']-d.iloc[1]['Lambda'])
                t = pd.DataFrame([[state,0.00,l0,0.00]], columns=['state','Lambda','DVDL', 'rms'])
                df = df.append(t)
            except:
                print('Error on the interpolation...most like  missing point')

        if d[d['Lambda'] == 1.0].empty:

            print("WARNING: Lambda 1 on {}, empty, data point will be interpolated".format(state))
            
            try:
                x = ((d.iloc[-2]['Lambda']-1.0)*d.iloc[-1]['DVDL'])+((1.0-d.iloc[-1]['Lambda'])*d.iloc[-2]['DVDL'])
                y = (d.iloc[-2]['Lambda']-d.iloc[-1]['Lambda'])
                t = pd.DataFrame([[state,1.00,x/y,0.00]], columns=['state', 'Lambda','DVDL', 'rms'])
                df = df.append(t)  
            except:                
                print('Error on the interpolation...most like  missing point')



    return df

def calc_ddg(df, STATES, output_file):

    delta_file = open(output_file, 'w')
    final = dict()
    for idx,state in enumerate(STATES):
        i = list()
        report = df[(df['state']==state)]
        integration = get_integration(report)
        i.append(integration)
        print('dG sum for {}:  \t\t {:.6f}'.format(state,sum(i)))
        print('dG sum for {}:  \t\t {:.6f}'.format(state,sum(i)), file=delta_file)
        print('-'*30)
        final[state] = sum(i)

    print('#'*30)
    # ΔGcomplex - ΔGligands
    ddg = final['complex']-final['ligand']
    print('final ddG binding: \t\t {:.6f}'.format( ddg))
    print('final ddG binding: \t\t {:.6f}'.format( ddg ), file= delta_file)

    delta_file.close()
    return ddg

def plot_integration(df, STATES, output_file):

    fig, ax = plt.subplots(nrows=2, ncols=1, figsize=(6,10))
    
    for idx,state in enumerate(STATES):


        report = df[(df['state']==state)]
        report.sort_values('Lambda', inplace=True)
        ax[idx].errorbar(report['Lambda'], report['DVDL'], yerr=report['rms'])
        ax[idx].set_title(state)

    fig.savefig(output_file)
    plt.close(fig)

def parse_TI(STATES, args):

    # Arguments
    # if is not none, change type
    if args.limit:
        LIMIT = int(args.limit)
        if LIMIT<args.skip:

            print('ERROR: Limit ({}) is smaller than the skip ({}) value '.format(LIMIT, args.skip))
            sys.exit(1)
    else:
        LIMIT = args.limit

    # step-lambda format# ex. complex-0.100
    if args.ignore:
        ignore_state, ignore_lambda = args.ignore.split('-')
    else:
        ignore_state = None
        ignore_lambda = None

    # Parse TI output
    preprocess_dvdl = list()
    for state in STATES:
        # setup ignore
        if ignore_state == state:
            IGNORE = ignore_lambda
        else:
            IGNORE = ''

        # get lamdba folders
        folders = glob.glob('./'+state+'/*')
        labels=[state]
        parse_DVDL_MD(folders, labels , preprocess_dvdl, skip=args.skip, limit=LIMIT, ignore=IGNORE, no_plot=args.no_plot)


    df = pd.DataFrame(preprocess_dvdl, columns=['state' ,'Lambda','DVDL', 'rms'])

    return df
 
def report(df, STATES):

    print('report:')
    for idx,state in enumerate(STATES):
        report = df[(df['state']==state)]
        report.sort_values('Lambda', inplace=True)
        print(' {} : {} Lamdbas'.format(state,  report['Lambda'].shape[0]))

        print(' Using : {} '.format( (report['Lambda'].unique())))

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
# skip multiple lambdas
# improve RMSD generation, may be use pjtraj to get values


def get_args():
    ## todo simplify mutations to X#X
    parser = argparse.ArgumentParser()
    parser.add_argument("--no-extrapolation", dest= "extrapol",default=True, action='store_false', help="by default if lambda 0 or 1 is missing the value is extrapolated. this argument disables this behaibour")
    parser.add_argument("--ignore", type=str, default=None, help="Ignore a lambda, format ->  step-lambda e.g.. complex-0.100 " )
    parser.add_argument("--skip", type=int, default=0, help='skip # of step to calculate the average on each lambda, 500 or 1ns recommended value, default = 0' )
    parser.add_argument("--limit",  default=None, help='limit on the  # of step to calculate the average dvdl for each lambda, default until the end'  )
    parser.add_argument("--no-plot", default=False, dest='no_plot', action='store_true', help='do not generate plots, only the ddg' )

    args = parser.parse_args()
    return args



if __name__ == '__main__':

    ## CONSTANTS
    STATES = ['complex', 'ligand']

    # get arguments
    args = get_args()

    # make folder
    create_folder('./analysis')
    if not args.no_plot:
        create_folder('./analysis/lambdas')
        create_folder('./analysis/rmsd')

    # Parse TI output
    df = parse_TI(STATES, args)

    # add lambda 0 and 1, lineal extrapolation
    if args.extrapol:
        df = extrapolation(df, STATES)

    # save lambdas
    df.to_csv('./analysis/dvdldata.csv', index=False)

    # if verbose
    # print report of lambdas
    report(df, STATES)

    # Calc ddg priint and write
    print('\n')
    calc_ddg(df, STATES, output_file='./analysis/delta.out')
    # plot integrations
    if not args.no_plot:
        plot_integration(df, STATES, output_file='./analysis/plot_integrations.png')

    # run it through pytraj
    # Plot RMSD
    for state in STATES:
        folders = glob.glob('./'+state+'/*')
        raw_H_data = defaultdict(list)

        for clambda in folders:
            if os.path.isfile(clambda+'/rmsd_bb.1.dat'):
                l = pd.read_csv(clambda+'/rmsd_bb.1.dat', delim_whitespace=True)
                s = os.path.basename(os.path.normpath(clambda))
                plot_rmsd(l,s, state)
