from __future__ import print_function
"""
.
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


from protocols import STD_INPUT, GEN_INPUT

AmberProtocols = GEN_INPUT

### All these protocols should be a class, TODO
# class AmberProtocols(object):
#     """docstring for ."""
#     def __init__(self, clambdas, masks, **args ):
#
#         self.clambas = clambdas
#         self.maks
def run_parmed(folder, resid, len_lig, system):
    '''Deduplication and mask
    '''
    # protein.parm7
    ligan_parm = """loadRestrt {system}_vdw_bonded.rst7
                    setOverwrite True
                    tiMerge :1-{} :{}-{} :{} :{}
                    outparm {system}_vdw_bonded.parm7 {system}_vdw_bonded.rst7
                    outpdb merged_{system}.pdb
                    quit""".format(len_lig, len_lig + 1, len_lig*2, resid, int(resid)+len_lig, system=system)
    print(ligan_parm, file=open(folder+'/{}.parmed'.format(system),'w') )

    process = subprocess.Popen('parmed {0}_vdw_bonded.parm7 -i {0}.parmed'.format(system).split(), cwd=folder, stdout=subprocess.PIPE)
    # wait for the process to terminate
    out = process.communicate()
    # both process return same mask

    masks = parse_masks_parmed(out)

    return masks

def parse_masks_parmed(output):

    mask = defaultdict(list)

    s = output[0].decode("utf-8")
    for line in s.split('\n'):
        if line.startswith('timask'):
            mask['timask'].append(line.strip())
        # if line.startswith('scmask'):
        #     mask['scmask'].append(line.strip())


    return mask

def compile_mask(status, masks, pos, l):

    #mask = defaultdict(list)
    mask = masks.copy()
    # mask assuming Ligand are first in the pdb
    mask['scmask'].append('''scmask1=':{}& !@CA,C,O,N,CB,H,HA', '''.format(pos))
    mask['scmask'].append('''scmask2=':{}& !@CA,C,O,N,CB,H,HA', '''.format(l+1))



    if status == 'vdw_bonded':
        # setup decharge
        mask['scmask'][1] = mask['scmask'][1]  + """ crgmask=':{}|:{}' ,""".format(pos,l+1)

    if status == 'decharge':

        mask['scmask'][1] = mask['scmask'][1]  + """ crgmask=':{}' ,""".format(l+1)

    if status == 'recharge':

        mask['scmask'][1] = mask['scmask'][1]  + """ crgmask=':{}' ,""".format(pos)

    return mask

def create_amber_scripts(folder, l, masks, increment, msteps=2000, hsteps=5000, psteps=40000):



#    for l in clambdas:
        for step, proto in AmberProtocols.items():
            # fix steps, using a class
            prot_input = proto.format(clambda=l, msteps=msteps, hsteps=hsteps, psteps=psteps,
                                timask1=masks['timask'][0],
                                timask2=masks['timask'][1],
                                scmask1=masks['scmask'][0],
                                scmask2=masks['scmask'][1],
                                increment=increment)
            path = '{}/{:.3f}/{}.in'.format(folder,  l, step)
            # MAIN+'/{}/{}/{:.3f}/ti.rst7'.format(s,state,w)
            output = open(path,'w')
            print(prot_input, file=output)
            output.close()



def mutate(res, resid, chain, mutatation, template, output):
    """Quick and diry mutation function. TODO change to biopython or something much robust.

    Parameters
    ----------

    Returns
    -------

    """
    residues = dict()
    backbone = ['C','O', 'N', 'CA', 'CB', 'H', 'HA']
    with open(template,'r') as input_file:
        for line in input_file:
            if line.startswith('ATOM') or line.startswith('HETATM'):
                atom_data = {'label':str(line[0:6].strip()),
                             'atom_num':int(line[6:11].strip()),
                             'atom_type':line[12:16].strip(),
                             'alternate':line[16:17].strip(),
                             'res_type':line[17:20].strip(),
                             'chain':line[21:22].strip(),
                             'res_num':int(line[22:26].strip()),
                             'insertion':line[26:27].strip(),
                             'x_coord':float(line[30:38].strip()),
                             'y_coord':float(line[38:46].strip()),
                             'z_coord':float(line[46:54].strip()),
                             'occupancy':float(line[54:60].strip()),
                             'tempFactor':float(line[60:66].strip()),
                             'element_symbol':line[76:78].strip(),
                             'charge':line[78:80].strip(),
                             'ASAS':'',
                             'vdw':'',
                             'hbond_type':'',
                             'coulomb':'',
                             'PE':''
                             }
                residues[atom_data['res_num']] = 1
            else:
                print(line.rstrip(), file=output)
                continue

            if atom_data['res_type'] == res  and atom_data['res_num'] == int(resid):
                if  atom_data['atom_type'] in backbone:
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



def merge_topologies(folder):

    tleap_input = """
    # load the AMBER force fields

    source leaprc.protein.ff14SB
    source leaprc.gaff
    source leaprc.water.tip3p
    loadAmberParams frcmod.ionsjc_tip3p



    # load the coordinates and create the systems
    # ligand = loadpdb $basedir/bnz.pdb
    m1 = loadpdb lig.pdb
    m2 = loadpdb mut_lig.pdb
    target = loadpdb trgt.pdb

    ligand = combine {m1 m2}
    complex = combine {m1 m2 target}

    # do not recenter coordinates
    set default nocenter on

    # create ligand in solution
    Addions ligand NA 0
    Addions ligand CL 0
    #Addions ligand NA 1
    #Addions ligand CL 1
    #solvateOct ligand TIP3PBOX 12.0
    solvateBox ligand TIP3PBOX 12.0

    savepdb ligand ligand_vdw_bonded.pdb
    saveamberparm ligand ligand_vdw_bonded.parm7 ligand_vdw_bonded.rst7

    # create complex in solution
    Addions complex NA 0
    Addions complex CL 0
    #Addions complex NA 1
    #Addions complex CL 1
    #solvateOct complex TIP3PBOX 12.0
    solvateBox complex TIP3PBOX 12.0

    savepdb complex complex_vdw_bonded.pdb
    saveamberparm complex complex_vdw_bonded.parm7 complex_vdw_bonded.rst7

    quit
    """

    print(tleap_input, file=open(folder+'/merge_tleap.inp','w') )

    process = subprocess.Popen('tleap -f merge_tleap.inp'.split(), cwd=folder)
    # wait for the process to terminate
    #out = process.communicate()
    #errcode = process.returncode
    status = process.wait()
    return

def create_topologies(folder):

    tleap_input = """
    # load the AMBER force fields
    source leaprc.protein.ff14SB
    source leaprc.gaff
    source leaprc.water.tip3p
    loadAmberParams frcmod.ionsjc_tip3p

    # coordinates for solvated ligands as created previously by MD
    lsolv = loadpdb ligand_solvated.pdb
    lwt = loadpdb ligand_wt.pdb
    lmut = loadpdb ligand_mut.pdb

    # coordinates for complex as created previously by MD
    csolv = loadpdb complex_solvated.pdb
    cwt = loadpdb complex_wt.pdb
    cwt = loadpdb complex_mut.pdb

    # decharge transformation
    decharge = combine {lwt lwt lsolv}
    setbox decharge vdw
    savepdb decharge ligand_decharge.pdb
    saveamberparm decharge ligand_decharge.parm7 ligand_decharge.rst7

    decharge = combine {cwt cwt csolv}
    setbox decharge vdw
    savepdb decharge complex_decharge.pdb
    saveamberparm decharge complex_decharge.parm7 complex_decharge.rst7

    # recharge transformation
    recharge = combine {lmut lmut lsolv}
    setbox recharge vdw
    savepdb recharge ligand_recharge.pdb
    saveamberparm recharge ligand_recharge.parm7 ligand_recharge.rst7

    recharge = combine {cwt cwt csolv}
    setbox recharge vdw
    savepdb recharge complex_recharge.pdb
    saveamberparm recharge complex_recharge.parm7 complex_recharge.rst7

    quit
    """

    print(tleap_input, file=open(folder+'/create_tleap.inp','w') )

    process = subprocess.Popen('tleap -f create_tleap.inp'.split(), cwd=folder)
    # wait for the process to terminate
    #out = process.communicate()
    #errcode = process.returncode
    status = process.wait()
    return


def strip_comp(folder, system, resid, len_lig):
    mut_start = len_lig+1
    mut_end = len_lig+1

    trajin = """parm {system}_vdw_bonded.parm7
                trajin {system}_vdw_bonded.rst7

                # remove the two ligands and keep the rest
                strip ":1-{mut_end}"
                outtraj {system}_solvated.pdb onlyframes 1

                # extract the first ligand
                unstrip
                strip ":{mut_start}-999999"
                outtraj {system}_wt.pdb onlyframes 1

                # extract the second ligand
                unstrip
                strip ":1-{len_lig},{mut_end}-999999"
                outtraj {system}_mut.pdb onlyframes 1

                """.format(system=system, len_lig = len_lig, mut_start= mut_start, mut_end = mut_end)

    print(trajin, file=open(folder+'/{0}_vdw_bonded.cppinp'.format(system),'w') )
    process = subprocess.Popen('cpptraj -i {0}_vdw_bonded.cppinp '.format(system).split(), cwd=folder, stdout=subprocess.PIPE)
        # wait for the process to terminate
    out = process.communicate()




def create_folders_tree(MAIN, LIG, COM, clambdas):


    mkdir(MAIN+'/'+LIG)
    mkdir(MAIN+'/'+COM)
    for state in ['decharge' ,'vdw_bonded' ,'recharge']:
        for window in clambdas:
            path =  '{}/{}/{}/{:.3f}'.format(MAIN,LIG, state,   window)
            mkdir(path)
            path =  '{}/{}/{}/{:.3f}'.format(MAIN,COM, state, window)
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
    ## todo simplify mutations to X#X
    parser = argparse.ArgumentParser()
    parser.add_argument("--ligand", type=str)
    parser.add_argument("--chain", type=str, default='')
    parser.add_argument("--mutation", type=str )
    parser.add_argument("--target", type=str )
    parser.add_argument("--increment", type=float, default=.1 )
    parser.add_argument("--msteps", type=int,default=4000 )
    parser.add_argument("--hsteps", type=int,default=5000 )
    parser.add_argument("--psteps", type=int,default=2500000 )
    #parser.set_defaults(nice=False)
    args = parser.parse_args()
    return args


if __name__ == '__main__':
    args = get_args()
    # generate lambdas
    res_num = args.mutation[3:-3]
    res_wt =  args.mutation[:3]
    res_mut = args.mutation[-3:]
    windows = np.arange(0,1+args.increment,args.increment)
    # Create Folders  & copy data
    MAIN = './MUTATION/' +  args.mutation
    # systems
    LIG = 'ligand'
    COM = 'complex'
    systems = [LIG, COM]


    create_folders_tree(MAIN, LIG, COM, windows)
    shutil.copy(args.ligand, MAIN+'/lig.pdb' )
    shutil.copy(args.target, MAIN+'/trgt.pdb' )
    # Simple MUTATION
    mut_ligand = open(MAIN+'/mut_lig.pdb','w')
    residues = mutate(res_wt, res_num, args.chain, res_mut, MAIN+'/lig.pdb', mut_ligand)
    mut_ligand.close()
    # merge tleap
    merge_topologies(MAIN)
    # strip water, and molecules
#    for s in systems:
#        strip_comp(MAIN, s, args.resid, len(residues))

    # Create topologies
    #create_topologies(MAIN)

    # run parmed
    for s in systems:
        parsed_masks = run_parmed(MAIN, res_num, len(residues), s)



    for s in systems:
        for state in ['decharge' ,'vdw_bonded' ,'recharge']:
            masks = compile_mask(state, parsed_masks, res_num, len(residues) )
            # masks control
            if len(masks) == 0:
                print('masks empty')
                raise ValueError
            if len(masks['scmask']) != 2 or len(masks['timask']) != 2:
                print('masks incomplete')
                print(masks)
                raise ValueError

        #create_amber_scripts(MAIN, s, windows, masks, steps=args.steps)
            for w in windows:
                shutil.copy(MAIN+'/{}_vdw_bonded.parm7'.format(s), MAIN+'/{}/{}/{:.3f}/ti.parm7'.format(s,state,w) )
                shutil.copy(MAIN+'/{}_vdw_bonded.rst7'.format(s), MAIN+'/{}/{}/{:.3f}/ti.rst7'.format(s,state,w) )
                #  create_amber_scripts(folder, l, masks, increment
                fol = MAIN+'/{}/{}/'.format(s,state)
                create_amber_scripts(fol, w, masks, increment=args.increment, msteps=args.msteps, hsteps=args.hsteps, psteps=args.psteps)


    # create restart script



    #
    #     s = output[0].decode("utf-8")
    #     for line in s.split('\n'):
    #         if line.startswith('timask'):
    #             mask['timask'].append(line.strip())
    #         if line.startswith('scmask'):
    #             mask['scmask'].append(line.strip())
    #
    #
    #     return mask
