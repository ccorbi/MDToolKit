import sys
import os
import glob


#def step_status

def get_steps(output_tis):
    last = 0

    for output in output_tis:
        buffer = 0
        with open(output) as fn:
            for line in fn:
                try:
                    if line.startswith(' NSTEP ='):
                
                        info = line.split()
                        buffer = int(info[2])
                except ValueError:
                    print('WARNING!, Parsing error on {}'.format(output))
                    pass

        last += buffer

    return last


def check_preti(lfolder):

    for pt_output in ['min.1.out', 'min.2.out','equil.3.out', 'equil.4.out']:
        if os.path.exists(lfolder+'/'+pt_output):
            pass

def print_status(lfolder):


    ti_outputs = glob.glob(lfolder+'/ti.*.out')
    ti_outputs.sort()
    total_steps = get_steps(ti_outputs)
    l = os.path.basename(lfolder)
    print(l, total_steps)

cwd = os.getcwd()

for status_folder in ['complex','ligand', 'folded', 'unfolded']:
    if os.path.isdir(cwd+'/'+status_folder):
        print(status_folder)
        lambdas = glob.glob(cwd+'/'+status_folder+'/*.*0')
        lambdas.sort()
        for l in lambdas:
            print_status(l)