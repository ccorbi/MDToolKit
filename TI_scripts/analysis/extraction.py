from __future__ import print_function
from __future__ import division
import math, sys
import os.path
import glob


class dVdL(object):
  '''A class that uses an online algorithm to compute mean and variance.'''

  def __init__(self, store_data = False):
    self.step = 0.0
    self.mean = 0.0
    self.M2 = 0.0

    self.store = store_data
    self.data = []


  def accumulate(self, x):
    '''Accumulate data points to compute mean and variance on-the-fly.'''

    self.step += 1

    delta = float(x) - self.mean

    self.mean += delta / float(self.step)
    self.M2 += delta * (float(x) - self.mean)

    if self.store:
      self.data.append(x)


  def get_variance(self):
    '''Convenience funtion to return variance.'''

    return self.M2 / (self.step - 1)


  def get_stat(self):
    '''Convenience funtion to return mean and standard deviation.'''

    return self.mean, math.sqrt(self.M2 / (self.step - 1))


#to use the script, choose the right fn number and the right filename prefix
#################################################################################################################
## loops over "fn" of output files .
## fn = file number
def read_TIMD_out(fileprefix,lbda):

    ### get average and rms values of DV/DL, TEMP, PRESS, Density parameters

    average_dvdl=list()
    rms_dvdl=list()
    for i in range(1, lbda+1):
        dv = dVdL()
        # collect reboots
        lambda_fileprefix = fileprefix.replace('prod','prod*' )
        lambda_fileprefix = lambda_fileprefix + str(i) + '.out'

        lambda_files = glob.glob(lambda_fileprefix)
        for l in lambda_files:
            with open(l,'r') as input_file:
                for line in input_file:

                    if 'R M S' in line:
                        break
                    if 'DV/DL' in line:
                        i = line.split()[-1]
                        dv.accumulate(i)


        a, s = dv.get_stat()
        average_dvdl.append(a)
        rms_dvdl.append(s)



    return average_dvdl,rms_dvdl
#################################################################################################################
def write_TIMD_report(ksteps,average_dvdl,rms_dvdl, lbda, fileout):


    fileout.write("Lambda  DVDL   rms   ksteps TEMP  rms    PRESS  rms   Density rms"+"\n")
    for i in range(1,lbda+1):
        fileout.write( "%5.3f\t%9.3f\t%6.3f" % (i/10.0, average_dvdl[i-1],rms_dvdl[i-1] ) )
        fileout.write( "\n" )
    return
#################################################################################################################
def main():

    #fileprefix="native-lig_prod_v0_l";
    fileprefix=sys.argv[1]
    lbda=9##number of files to be read
    ksteps=500000

    average_dvdl, rms_dvdl = read_TIMD_out( fileprefix, lbda )

    file_report = "report." + fileprefix
    fileout = open( file_report, 'write' )
    write_TIMD_report(ksteps,average_dvdl,rms_dvdl, lbda, fileout)
    fileout.close()
#################################################################################################################
#################################################################################################################
main()
