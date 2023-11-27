import lhsmdu
import numpy as np
import pandas as pd
from math import sqrt
  
lhsmdu.setRandomSeed(None)

delta_diago = 0.10
delta_shear = 0.10
min_diago = 1. - 0.5*delta_diago
min_shear =    - 0.5*delta_shear
eeq_threshold = 0.1
I33 = np.identity(3)

def CG(F11,F22,F33,F12,F13,F23):

  F = np.array([[F11, F12, F13],
		[  0, F22, F23],
		[  0,   0, F33]])

  return np.dot(np.transpose(F), F)

def E(F11,F22,F33,F12,F13,F23):
  return 0.5 * (CG(F11,F22,F33,F12,F13,F23) - I33)

def Edev(F11,F22,F33,F12,F13,F23):
  e = E(F11,F22,F33,F12,F13,F23)
  tre = np.trace(e)
  return (e - I33 * tre / 3.)

def Eeq(F11,F22,F33,F12,F13,F23):
  deve = Edev(F11,F22,F33,F12,F13,F23)
  return sqrt(2 * np.sum(deve*deve) / 3.)

npts=500
nvars=6

liste = pd.DataFrame((lhsmdu.sample(nvars,npts)).T)
liste.columns = ['F11', 'F22', 'F33', 'F12', 'F13', 'F23']

liste['F11'] = min_diago + delta_diago*liste['F11']
liste['F22'] = min_diago + delta_diago*liste['F22']
liste['F33'] = min_diago + delta_diago*liste['F33']

liste['F12'] = min_shear + delta_shear*liste['F12']
liste['F13'] = min_shear + delta_shear*liste['F13']
liste['F23'] = min_shear + delta_shear*liste['F23']
    
for index, row in liste.iterrows():
  eeq = Eeq(row['F11'],row['F22'],row['F33'],row['F12'],row['F13'],row['F23'])
  if eeq > eeq_threshold:
    liste.drop(index, inplace=True)
liste_np = liste.to_numpy()

nptsmax=100
for j in range(nptsmax):
  ofile = open('dir.deformation.all/deformation.%d.mod' %(j), 'w')
  ofile.write('change_box     all triclinic\n')
  ofile.write('variable       F11 equal %5.4f\n' %(liste_np[j,0]))
  ofile.write('variable       F22 equal %5.4f\n' %(liste_np[j,1]))
  ofile.write('variable       F33 equal %5.4f\n' %(liste_np[j,2]))
  ofile.write('variable       F12 equal %5.4f\n' %(liste_np[j,3]))
  ofile.write('variable       F13 equal %5.4f\n' %(liste_np[j,4]))
  ofile.write('variable       F23 equal %5.4f\n' %(liste_np[j,5]))    
  ofile.write('variable       lxnew equal ${F11}*lx\n')
  ofile.write('variable       xynew equal ${F12}*ly\n')
  ofile.write('variable       lynew equal ${F22}*ly\n')    
  ofile.write('variable       xznew equal ${F13}*lz\n')
  ofile.write('variable       yznew equal ${F23}*lz\n')    
  ofile.write('variable       lznew equal ${F33}*lz\n')
  ofile.write('change_box     all x final 0.0 ${lxnew} y final 0.0 ${lynew} z final 0.0 ${lznew} xy final ${xynew} xz final ${xznew} yz final ${yznew} remap units box\n')
  
  # |F11  F12  F13| |xmax-xmin  0          0        | |F11*(xmax-xmin)  F12*(ymax-ymin) F13*(zmax-zmin)|
  # |0    F22  F23|.|0          ymax-ymin  0        |=|0                F22*(ymax-ymin) F23*(zmax-zmin)|
  # |0    0    F33| |0          0          zmax-zmin| |0                0               F33*(zmax-zmin)|
  F=np.zeros((3,3))
  F[0,0]=liste_np[j,0]
  F[1,1]=liste_np[j,1]
  F[2,2]=liste_np[j,2]
  F[0,1]=liste_np[j,3]
  F[0,2]=liste_np[j,4]
  F[1,2]=liste_np[j,5]
    
  ofile.close()
  
