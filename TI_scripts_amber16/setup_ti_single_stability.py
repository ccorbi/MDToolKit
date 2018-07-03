#!/usr/bin/env python

from __future__ import print_function
"""
.based on https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2872215/
"""
import sys
import os
import argparse
import subprocess
import glob
import re
import shutil
from collections import defaultdict
import numpy as np

from protocols import STD_INPUT, GEN_INPUT, STB_INPUT

AmberProtocols = STB_INPUT

# All these protocols should be a class, TODO
# class AmberProtocols(object):
#     """docstring for ."""
#     def __init__(self, clambdas, masks, **args ):
#
#         self.clambas = clambdas
#         self.maks


def compile_mask(pos, l, mask, ignore):

    mask['scmask'].append('''scmask1=':{} {}', '''.format(pos, ignore))  # CB
    mask['scmask'].append('''scmask2=':{} {}', '''.format(l + 1, ignore))  # CB

    return mask


def create_amber_scripts(folder, system, clambdas, masks, increment, psteps=40000):

    for l in clambdas:
        for step, proto in AmberProtocols.items():
            # fix steps, using a class
            prot_input = proto.format(clambda=l, psteps=psteps,
                                      timask1=masks['timask'][0],
                                      timask2=masks['timask'][1],
                                      scmask1=masks['scmask'][0],
                                      scmask2=masks['scmask'][1],
                                      increment=increment)
            path = '{}/{}/{:.3f}/{}.in'.format(folder, system, l, step)
            output = open(path, 'w')
            print(prot_input, file=output)
            output.close()


def mutate(res, resid, chain, mutatation, template, output, fix_atoms):
    """Quick and diry mutation function. TODO change to biopython or something much robust.

    Parameters
    ----------

    Returns
    -------

    """
    residues = dict()

    with open(template, 'r') as input_file:
        for line in input_file:
            if line.startswith('ATOM') or line.startswith('HETATM'):
                atom_data = {'label': str(line[0:6].strip()),
                             'atom_num': int(line[6:11].strip()),
                             'atom_type': line[12:16].strip(),
                             'alternate': line[16:17].strip(),
                             'res_type': line[17:20].strip(),
                             'chain': line[21:22].strip(),
                             'res_num': int(line[22:26].strip()),
                             'insertion': line[26:27].strip(),
                             'x_coord': float(line[30:38].strip()),
                             'y_coord': float(line[38:46].strip()),
                             'z_coord': float(line[46:54].strip()),
                             'occupancy': float(line[54:60].strip()),
                             'tempFactor': float(line[60:66].strip()),
                             'element_symbol': line[76:78].strip(),
                             'charge': line[78:80].strip(),
                             'ASAS': '',
                             'vdw': '',
                             'hbond_type': '',
                             'coulomb': '',
                             'PE': ''
                             }
                residues[atom_data['res_num']] = 1
            else:
                print(line.rstrip(), file=output)
                continue

            if atom_data['res_type'] == res and atom_data['res_num'] == resid:
                if atom_data['atom_type'] in fix_atoms:
                    line = line.replace(res, mutatation)
                    # save
                    print(line.rstrip(), file=output)
                else:
                    continue
            else:
                # remove hidrogens
                # if not re.match(r'^H', atom_data['atom_type']):
                print(line.rstrip(), file=output)
    return residues



def flank_mutate(resid, chain, mutatation, template, output, fix_atoms):
    """Quick and diry mutation function. TODO change to biopython or something much robust.

    Parameters
    ----------

    Returns
    -------

    """
    residues = dict()

    with open(template, 'r') as input_file:
        for line in input_file:
            if line.startswith('ATOM') or line.startswith('HETATM'):
                atom_data = {'label': str(line[0:6].strip()),
                             'atom_num': int(line[6:11].strip()),
                             'atom_type': line[12:16].strip(),
                             'alternate': line[16:17].strip(),
                             'res_type': line[17:20].strip(),
                             'chain': line[21:22].strip(),
                             'res_num': int(line[22:26].strip()),
                             'insertion': line[26:27].strip(),
                             'x_coord': float(line[30:38].strip()),
                             'y_coord': float(line[38:46].strip()),
                             'z_coord': float(line[46:54].strip()),
                             'occupancy': float(line[54:60].strip()),
                             'tempFactor': float(line[60:66].strip()),
                             'element_symbol': line[76:78].strip(),
                             'charge': line[78:80].strip(),
                             'ASAS': '',
                             'vdw': '',
                             'hbond_type': '',
                             'coulomb': '',
                             'PE': ''
                             }
                residues[atom_data['res_num']] = 1
            else:
                print(line.rstrip(), file=output)
                continue

            if atom_data['res_num'] == resid:
                if atom_data['atom_type'] in fix_atoms:
                    line = line.replace(atom_data['res_type'], mutatation)
                    # save
                    print(line.rstrip(), file=output)
                else:
                    continue
            else:
                # remove hidrogens
                #if not re.match(r'^H', atom_data['atom_type']):
                print(line.rstrip(), file=output)
    return residues


def merge_topologies(folder):

    tleap_input = """
    # load the AMBER force fields
    source leaprc.protein.ff14SB
    source leaprc.gaff
    source leaprc.water.tip3p
    loadAmberParams frcmod.ionsjc_tip3p

    # load the coordinates and create the systems
    m1 = loadpdb unfolded.pdb
    m2 = loadpdb mut_unfolded.pdb
    target = loadpdb trgt.pdb
    target_mut = loadpdb  mut_trgt.pdb

    unfolded_state = combine {m1 m2}
    folded_state = combine {target target_mut}

    # do not recenter coordinates
    set default nocenter on

    # create unfolded_state in solution
    Addions unfolded_state NA 0
    Addions unfolded_state CL 0

    solvateBox unfolded_state TIP3PBOX 12.0
    savepdb unfolded_state unfolded.pdb
    saveamberparm unfolded_state unfolded.parm7 unfolded.rst7

    # create folded_state in solution
    Addions folded_state NA 0
    Addions folded_state CL 0

    solvateBox folded_state TIP3PBOX 12.0
    savepdb folded_state folded.pdb
    saveamberparm folded_state folded.parm7 folded.rst7

    quit
    """

    print(tleap_input, file=open(folder + '/merge_tleap.inp', 'w'))

    process = subprocess.Popen('tleap -f merge_tleap.inp'.split(), cwd=folder)
    # wait for the process to terminate
    #out = process.communicate()
    #errcode = process.returncode
    status = process.wait()
    return


def run_parmed(folder, resid, len_lig, system):
    '''Deduplication and mask
    '''
    # protein.parm7
    ligan_parm = """loadRestrt {system}.rst7
                    setOverwrite True
                    tiMerge :1-{} :{}-{} :{} :{}
                    outparm merged_{system}.parm7 merged_{system}.rst7
                    outpdb merged_{system}.pdb
                    quit""".format(len_lig, len_lig + 1, len_lig * 2, resid, int(resid) + len_lig, system=system)
    print(ligan_parm, file=open(folder + '/{}.parmed'.format(system), 'w'))

    process = subprocess.Popen(
        'parmed {0}.parm7 -i {0}.parmed'.format(system).split(), cwd=folder, stdout=subprocess.PIPE)
    # wait for the process to terminate
    out = process.communicate()
    # both process return same mask
    masks = parse_masks_parmed(out)
    # print(out)
    return masks


def parse_masks_parmed(output):

    mask = defaultdict(list)

    s = output[0].decode("utf-8")
    for line in s.split('\n'):
        if line.startswith('timask'):
            mask['timask'].append(line.strip())
        # if line.startswith('scmask'):
            # mask['scmask'].append(line.strip())

    return mask


def check_masks(masks):

    if len(masks) <2:
        print('Error masking')
        print(masks)
        sys.exit(1)
    try:
        for m in ['timask', 'scmask']:
            if len(masks[m]) <2:
                print('Error masking')
                print(masks)
                sys.exit(1)
    except:
        print('Error masking')
        print(masks)
        sys.exit(1)
    
    return

def create_folders_tree(MAIN, LIG, COM, clambdas):

    mkdir(MAIN + '/' + LIG)
    mkdir(MAIN + '/' + COM)
    for window in clambdas:
        path = '{}/{}/{:.3f}'.format(MAIN, LIG, window)
        mkdir(path)
        path = '{}/{}/{:.3f}'.format(MAIN, COM, window)
        mkdir(path)

    return


def mkdir(folder):

    try:
        os.makedirs(folder)
    except OSError as exc:
        # folder exist
        pass

    return


def get_args():
    # todo simplify mutations to X#X
    parser = argparse.ArgumentParser()
    parser.add_argument("--ligand", type=str)
    parser.add_argument("--chain", type=str, default='')
    parser.add_argument("--mutation", type=str)
    parser.add_argument("--target", type=str)
    parser.add_argument("--increment", type=float, default=.1)
    parser.add_argument("--scmask-ignore", type=str,
                        default='CA,C,O,N,HA,H1,H2,H3,H ', dest='scmask_ignore')
    parser.add_argument("--decouple-mask", 
                        default=False, dest='decouple_mask', action="store_true")
    parser.add_argument("--psteps", type=int, default=2500000)
    parser.add_argument("--ouput-folder", type=str,
                        default='MUTATION', dest='output_folder')
    # parser.set_defaults(nice=False)
    args = parser.parse_args()
    return args


if __name__ == '__main__':

    # load arguments
    args = get_args()
    # Parse mutation
    res_num = args.mutation[3:-3]
    res_wt = args.mutation[:3]
    res_mut = args.mutation[-3:]
    # generate lambdas
    windows = np.arange(0, 1 + args.increment, args.increment)
    # Create Folders  & copy data
    MAIN = './' + args.output_folder + '/' + args.mutation
    # systems Setup
    UNF = 'unfolded'
    COM = 'folded'
    systems = [UNF, COM]
    # create folders
    create_folders_tree(MAIN, UNF, COM, windows)
    #shutil.copy(args.target, MAIN + '/unfolded.pdb')
    shutil.copy(args.target, MAIN + '/trgt.pdb')

    # Parse&Setup masking
    fix_atoms = args.scmask_ignore.split(',')

    if 'PRO' in args.mutation:
        if 'H' in fix_atoms:
            print('WARNING remove H from the secondary mask')
            fix_atoms.remove('H')
            args.scmask_ignore = args.scmask_ignore.replace(',H ','')

    if args.decouple_mask:
        mask_ignore = ''
    else:
        mask_ignore = '& !@' + args.scmask_ignore

    # Simple MUTATION
    # add flanking GLYCINES
    mut_ligand = open(MAIN + '/unfolded_temp.pdb', 'w')
    residues = flank_mutate( int(res_num)-1, args.chain, 'GLY',
                      MAIN + '/trgt.pdb', mut_ligand, ['CA','C','O','N'])
    mut_ligand.close()

    mut_ligand = open(MAIN + '/unfolded.pdb', 'w')
    residues = flank_mutate( int(res_num)+1, args.chain, 'GLY',
                      MAIN + '/unfolded_temp.pdb', mut_ligand, ['CA','C','O','N'])
    mut_ligand.close()


    # create mutate unfolded
    mut_ligand = open(MAIN + '/mut_unfolded.pdb', 'w')
    residues = mutate(res_wt, int(res_num), args.chain, res_mut,
                      MAIN + '/unfolded.pdb', mut_ligand, fix_atoms)
    mut_ligand.close()

    # create mutate folded
    mut_ligand = open(MAIN + '/mut_trgt.pdb', 'w')
    residues = mutate(res_wt, int(res_num), args.chain, res_mut,
                      MAIN + '/trgt.pdb', mut_ligand, fix_atoms)
    mut_ligand.close()

    # merge tleap
    merge_topologies(MAIN)

    # run parmed 
    for s in systems:
        masks = run_parmed(MAIN, res_num, len(residues), s)
        print('output run_parmed')
        print(masks)
    # build masking
        masks = compile_mask(res_num, len(residues), masks, mask_ignore)
        print('afet')
        print(masks)
        check_masks(masks)


    # create amber scripts
    #for s in systems:
        create_amber_scripts(MAIN, s, windows, masks,
                             increment=args.increment, psteps=args.psteps)
        for w in windows:
            shutil.copy(MAIN + '/merged_{}.parm7'.format(s),
                        MAIN + '/{}/{:.3f}/ti.parm7'.format(s, w))
            shutil.copy(MAIN + '/merged_{}.rst7'.format(s),
                        MAIN + '/{}/{:.3f}/ti.rst7'.format(s, w))
