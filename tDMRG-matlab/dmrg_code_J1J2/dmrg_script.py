import os
import shutil
import numpy

#os.system("pdftk *.pdf.pdf shuffle output collate.pdf compress")
#os.system("convert collate.pdf -compress Zip collate_zip.pdf")

import subprocess

U=0
J=1
mu=0
C=10
N=10
Nph=1
NT=Nph
mBt=2
V0=1e-2
Nop=7

dt=0.1
nsweep=6
ntsweep=0
isweep=0
tisweep=0

ctrlsweep=4
convergencE=1e-12

def create_input():

  filein=open('initial.m','w',0)
  # filein.flush()
  filein.write('%20s\n\n' % ( 'function initial()' ) )
  filein.write('%20s\n\n' % ( 'global U J mu  V0 C N Nph NT' ) )
  filein.write('%20s\n\n' % ( 'global Nop ms Sops SopL SopR SEs' ) )
  filein.write('%20s\n\n' % ( 'global datafile isweep tisweep ntsweep ti' ) )
  filein.write('%2s %5s %1s\n' % ( 'U=', U, ';') )
  filein.write('%2s %5s %1s\n' % ( 'J=', J, ';') )
  filein.write('%2s %5s %1s\n' % ( 'mu=', mu, ';') )
  filein.write('%2s %5s %1s\n' % ( 'C=', C, ';') )
  filein.write('%2s %5s %1s\n' % ( 'N=', N, ';') )
  filein.write('%2s %5s %1s\n' % ( 'Nph=', Nph, ';') )
  filein.write('%2s %5s %1s\n' % ( 'mBt=', mBt, ';') )
  filein.write('%2s %5s %1s\n' % ( 'V0=', V0, ';') )
  filein.write('%2s %5s %1s\n' % ( 'Nop=', Nop, ';') )
  
  filein.write('%2s %5s %1s\n' % ( 'dt=', dt, ';') )
  filein.write('%2s %5s %1s\n' % ( 'nsweep=', nsweep, ';') )
  filein.write('%2s %5s %1s\n' % ( 'ntsweep=', ntsweep, ';') )
  filein.write('%2s %5s %1s\n' % ( 'isweep=', isweep, ';') )
  filein.write('%2s %5s %1s\n' % ( 'tisweep=', tisweep, ';') )
  filein.write('%2s %5s %1s\n' % ( 'ctrlsweep=', ctrlsweep, ';') )
  filein.write('%2s %5s %1s\n' % ( 'convergencE=', convergencE, ';') )

for J in numpy.arange(1,3,1):

  print J

  create_input()

  subprocess.call('octave dmrg_main.m', shell=True)
