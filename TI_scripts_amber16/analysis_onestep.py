#!/usr/bin/env python
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
import pytraj as pt
warnings.filterwarnings("ignore")




def get_rms(traj_path):

    traj = load_traj(traj_path)

    rmsd = pt.rmsd(traj, ref=0, mask='@CA,C,N')
    df_rmsd = pd.DataFrame(rmsd, columns=['RMSD'])
    
    rmsf = pt.rmsf(traj)
    df_rmsf = pd.DataFrame(rmsf, columns=['ATOM','RMSF'])

    return df_rmsd, df_rmsf


def load_traj(traj_path):
    ########## TAKE A LOOK HERE THIS WIERD!!!!
    # clamba do no defined, and in a notebook iterload failded to to this
    traj = pt.iterload(clambda+'/ti.*.nc', top=clambda+'/ti.parm7', frame_slice=(0,-1),)
    traj.autoimage()
    traj.superpose()

    return traj


def plot_rmsd(rmsd, clambda, state):

    fig, ax =  plt.subplots()

    rmsd['acc_RMSD'] = rmsd.expanding(2).mean()
    # free peptides trend to show higher  rmsd
    if state == 'ligand':
        ax.set_ylim(0,6)    
    else:
        ax.set_ylim(0,3)
    ax.plot(rmsd['RMSD'])
    ax.plot(range(rmsd.shape[0]), rmsd['acc_RMSD'], "r")
    #y_av = movingaverage(dvdl, 200)
    #ax.plot(range(len(dvdl)), y_av,"r")
    ax.set_title(clambda+ ' '+ state)
    fig.savefig('./analysis/rms/rmsd_{}_{}.png'.format(state,clambda ))
    plt.close(fig)

    return 



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
        if islambdafolder(l) and l not in ignore :
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



def read_ignore_lambdas(args):

    ignore = defaultdict(list)
    for ignore_point in args:
        state, lambdas = ignore_point.split('-')
        ignore[state].append(lambdas)

    return ignore


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
        #for ign_point in args.ignore
        ignore_points = read_ignore_lambdas(args.ignore)
    else:
        ignore_points = dict()

    # Parse TI output
    preprocess_dvdl = list()
    for state in STATES:
        # setup ignore
        if state in ignore_points:
            IGNORE = ignore_points[state]
        else:
            IGNORE = list()

        # get lamdba folders
        folders = glob.glob('./'+state+'/*')
        labels=[state]
        # Parse TI output
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
    parser.add_argument("--ignore", nargs='*', help="Ignore a lambda, format ->  step-lambda e.g.. complex-0.100 " )
    parser.add_argument("--skip", type=int, default=0, help='skip # of step to calculate the average on each lambda, 500 or 1ns recommended value, default = 0' )
    parser.add_argument("--limit",  default=None, help='limit on the  # of step to calculate the average dvdl for each lambda, default until the end'  )
    parser.add_argument("--no-plots", default=False, dest='no_plot', action='store_true', help='do not generate plots, only the ddg' )
    parser.add_argument("--no-rms", default=False, dest='no_rms', action='store_true', help='do not generate RMSD  and RMSF data' )

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
        create_folder('./analysis/rms')
        create_folder('./analysis/rms/data')

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
    if not args.no_rms:
        for state in STATES:
            folders = glob.glob('./'+state+'/*')
        #     raw_H_data = defaultdict(list)

            for clambda in folders:

                #if os.path.isfile(clambda+'/rmsd_bb.1.dat'):
                rmsd, rmsf = get_rms(clambda)
                s = os.path.basename(os.path.normpath(clambda))

                rmsd.to_csv('./analysis/rms/data/rmsd_{}_{}.csv'.format(s,state))
                rmsf.to_csv('./analysis/rms/data/rmsf_{}_{}.csv'.format(s,state))
                if not args.no_plot:
                    plot_rmsd(rmsd ,s, state)