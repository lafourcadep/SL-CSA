import sys
import numpy as np
import pandas as pd
import warnings
warnings.filterwarnings('ignore')
import pickle

def compute_ncomps(twojmax): # ok
    n_idx=0
    for j1 in range(0,twojmax+1):
        for j2 in range(0,j1+1):
            for j in range(int(abs(j1-j2)),min(twojmax,j1+j2)+1,2):
                if((j+j1+j2)%2 == 1):
                    continue
                if(j<j1):
                    continue
                n_idx+=1
    nidx=n_idx
    return nidx

specie=sys.argv[1]
cs=sys.argv[2]
nnn=int(sys.argv[3])
twojmax=int(sys.argv[4])
cstype=int(sys.argv[5])
ncomps_bs=compute_ncomps(twojmax)
bscomp_list = ['b%d'%i for i in range(1,ncomps_bs+1)]
header_list = ['CS'] + bscomp_list

ntfs=100
df_cs = pd.DataFrame(columns=header_list)
for k in range(ntfs):
    dump_file = 'bnnn_full_dumps_database/dump_%s_%s.%d.ramp'  %(specie,cs,k)
    data = np.loadtxt(dump_file, skiprows=9)
    df_loc = pd.DataFrame(data[:,4:] , columns=bscomp_list)
    cs_vector = [cstype for l in range(len(data))]
    df_loc.insert(0, "CS", cs_vector, True)
    df_cs = pd.concat([df_cs,df_loc])
stats_numeric = df_cs['CS'].describe()
print(stats_numeric)
df_cs.to_pickle('bnnn_full_pickle_database/database_bnnn_%s_%s_%d.pickle' %(specie,cs,cstype))
