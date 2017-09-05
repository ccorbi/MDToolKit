import math, sys
import os.path
from math import sqrt

#to use the script, choose the right filename prefix
#################################################################################################################
## loops over "fn" of output files .
## fn = file number 
def read_dvdl_out(filename):

    dvdl=[];rms_dvdl=[];
    file_op = open(filename, 'r');

    lines  =  file_op.readlines();
    counter = 0;

    for line in lines:
        linesplit = line.split() #split on white space
        if (counter > 0 and len(linesplit) > 0 ):
            dvdl.append(float(linesplit[1]));
            rms_dvdl.append(float(linesplit[2]));
        counter = counter + 1;
    file_op.close();
    return dvdl,rms_dvdl;
#################################################################################################################
def write_dvdl_prep( dvdl,rms_dvdl,file):

    length = len(dvdl);
    for i in range(1,length+1):
        file.write( "%5.3f\t%9.3f\t%6.3f" %(i/10.0,dvdl[i-1],rms_dvdl[i-1]) );
        file.write( "\n" );
    return;
#################################################################################################################
def main():

########calculate sub_dvdl and sub_rms_dvdl for different steps. 
    sub_dvdl=[];sub_rms_dvdl=[];

    for i in range(1,4):
        filename = "lig_step"+str(i)+".report";    
        l_dvdl,l_rms_dvdl=read_dvdl_out(filename);
        filename = "comp_step"+str(i)+".report";
        c_dvdl,c_rms_dvdl=read_dvdl_out(filename);
        for j in range(0,9):
            sub_dvdl.append( c_dvdl[j]-l_dvdl[j] );
            sub_rms_dvdl.append( sqrt( c_rms_dvdl[j]**2 + l_rms_dvdl[j]**2) );

        file_integ = "integ_prep_step"+str(i);
        file = open( file_integ, 'write' );
        write_dvdl_prep( sub_dvdl,sub_rms_dvdl,file );
        sub_dvdl = []; sub_rms_dvdl = []; 
        file.close()
#################################################################################################################
main()
