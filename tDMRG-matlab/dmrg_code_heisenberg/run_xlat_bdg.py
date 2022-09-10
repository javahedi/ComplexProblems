import os
import shutil
import numpy

#os.system("pdftk *.pdf.pdf shuffle output collate.pdf compress")
#os.system("convert collate.pdf -compress Zip collate_zip.pdf")

import subprocess

Ndim=2 
Nneigh=6
LL=[8, 8]
T=1.0e-2
U0=4.00             
Hartree='.false.'
OpenBC='.false.'
itmax=100
amix=1.0
eps=1.e-3
eps1=1.e-3
mu=-1.5 
HH=[0.4, 0.0, 0.0]
Ntot=100
Pol=[0.0, 0.0, 0.0]
alpha=[0.0, 0.0]


def create_input():

  filein=open('lattice_bdg.in','w',0)
  # filein.flush()

  filein.write('%4s %4s %20s\n' % ( Ndim, Nneigh                , '     ! Ndim, Nneigh') )
  filein.write('%8s     %20s\n' % ( ' '.join(map(str,LL))       , '     ! LL') )
  filein.write('%4s     %20s\n' % ( T                           , '     ! Temperature') )
  filein.write('%4s     %20s\n' % ( U0                          , '     ! U0') )
  filein.write('%4s     %20s\n' % ( Hartree                     , '     ! Hartree') )
  filein.write('%4s     %20s\n' % ( OpenBC                      , '     ! OpenBC ') )
  filein.write('%4s     %20s\n' % ( itmax                       , '     ! itmax') )
  filein.write('%4s     %20s\n' % ( amix                        , '     ! amix') )
  filein.write('%4s %4s %20s\n' % ( eps, eps1                   , '     ! eps, eps1') )
  filein.write('%4s %9s %20s\n' % ( mu, ' '.join(map(str,HH))   , '     ! mu, HH') )
  filein.write('%4s %9s %20s\n' % ( Ntot, ' '.join(map(str,Pol)), '     ! Ntot, Pol') )
  filein.write('%8s     %20s\n' % ( ' '.join(map(str,alpha))    , '     ! alpha') )


fileout=open('del_vs_U0.dat','w',0)

for U0 in numpy.arange(2,2.2,.2):

  create_input()

  subprocess.call('rm -f del_in.dat', shell=True)

  subprocess.call('time ./xlatbdg', shell=True)

  subprocess.call('cp bdg_soc2d.dat data/bdg_soc2d_U0_'+str(U0)+'.dat', shell=True)
  shutil.copy2('bdg_soc2d.dat', 'data/bdg_2d_U0_'+str(U0)+'.dat')

  for line in open('bdg_soc2d.dat'):
  	if "delmax" in line:
		fileout.write('%8s %8s \n' % (U0, line.split()[3]) )
