import numpy as np
import pandas as pd
import warnings
warnings.filterwarnings('ignore')
import pickle
import matplotlib as mpl
import matplotlib.pyplot as plt
from math import sqrt
from sklearn.metrics import classification_report, accuracy_score, f1_score
from sklearn.model_selection import train_test_split
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis as LDA
from sklearn.linear_model import LogisticRegression 
from sklearn.neighbors import KNeighborsClassifier
from sklearn.utils import shuffle
import sys
import os
import seaborn as sns
import joblib
import tkinter
mpl.use('TkAgg')
mpl.rc('text')#, usetex=True)
mpl.rcParams['ytick.labelsize']=30
mpl.rcParams['xtick.labelsize']=15
mpl.rcParams['axes.labelsize']=30
mpl.rcParams['legend.fontsize']=15
mpl.rcParams['axes.titlesize']=30
mpl.rcParams['figure.titlesize']=40


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

def class_mean_cov_full(X_train_lda, data_output, ncomps, path): #, extension):
  print('Mahalanobis distance calculation, inside each class')
  X = pd.concat([pd.DataFrame(X_train_lda), pd.DataFrame(data_output)], axis=1).dropna().sample(10000)
#  X.columns = ['1', '2', '3','CS'] 
  databs = X.values
  
  typessss=['bcc', 'fcc', 'hcp', 'dia']

  for i in range(len(typessss)):
    meanloc=np.zeros(ncomps)
    covloc=np.zeros((ncomps,ncomps))
    dataloc = databs[databs[:,ncomps] == 1.0*i][:,:ncomps].astype(float)

    for j in range(ncomps):
      meanloc[j] = np.mean(dataloc[:,j])                    
    ofile = open('%s/mean_database_%s.dat' %(path, typessss[i]), 'w')
    for j in range(ncomps):
      ofile.write('%5.8f\n' %(meanloc[j]))
    ofile.close()

    covloc = np.cov(dataloc)    
    ofile = open('%s/cov_database_%s.dat' %(path, typessss[i]), 'w')
    for j in range(ncomps):
      for k in range(ncomps):        
        ofile.write('%5.8f ' %(covloc[j,k]))
      ofile.write('\n')
    ofile.write('\n')
    ofile.close()

twojmax=int(sys.argv[1])
ncomps_bs=compute_ncomps(twojmax)
ncstypes=int(sys.argv[2])
cstypes = [i for i in range(ncstypes)]
crystal_structures = [sys.argv[3], sys.argv[5], sys.argv[7], sys.argv[ 9]]
crystal_types      = [sys.argv[4], sys.argv[6], sys.argv[8], sys.argv[10]]

ntfs=100
frames=[]
for i in range(ncstypes):
    pickle_name='bnnn_full_pickle_database/database_bnnn_%s_%s_%s.pickle' %(crystal_types[i], crystal_structures[i], cstypes[i])
    df_cs = pd.read_pickle(pickle_name)
    frames.append(df_cs)
df_full=pd.concat(frames).reset_index(drop=True)
df_train, df_test = train_test_split(df_full, test_size=0.2, shuffle=True)
stats_train = df_train['CS'].describe()
stats_test = df_test['CS'].describe()

# Training database
train_features = df_train.loc[:, ["b%d"%i for i in range(1,ncomps_bs+1)]]
train_labels   = df_train.loc[:, ['CS']]
train_features_data = train_features.values.astype(float)
train_labels_data   = train_labels.values.astype(int)

# Testing database
test_features = df_test.loc[:, ["b%d"%i for i in range(1,ncomps_bs+1)]]
test_labels   = df_test.loc[:, ['CS']]
test_features_data = test_features.values.astype(float)
test_labels_data   = test_labels.values.astype(int)

path='dir.slcsa'
isExist = os.path.exists(path)
if not isExist:
  os.makedirs(path)

##################################################################
# Linear Discriminant Analysis Start
print('Performing Linear Discriminant Analysis...')

lda_model = LDA(n_components=ncstypes-1)

train_features_data_reduced = lda_model.fit_transform(train_features_data, train_labels_data)
df_train_reduced = pd.concat([pd.DataFrame(train_features_data_reduced), pd.DataFrame(train_labels_data)], axis=1).dropna()
df_train_reduced.columns = ['C%d'%(i+1) for i in range(ncstypes-1)] + ['CS'] 

test_features_data_reduced  = lda_model.transform(test_features_data)
df_test_reduced = pd.concat([pd.DataFrame(test_features_data_reduced), pd.DataFrame(test_labels_data)], axis =1).dropna()
df_test_reduced.columns = ['C%d'%(i+1) for i in range(ncstypes-1)] + ['CS'] 

df_full_reduced = pd.concat([df_train_reduced,df_test_reduced])
stats_numeric = df_full_reduced['CS'].describe()

path='bnnn_reduced_pickle_database'
isExist = os.path.exists(path)
if not isExist:
  os.makedirs(path)

for i in range(ncstypes):
  df_cs = df_full_reduced.loc[df_full_reduced['CS'] == i]
  df_cs.to_pickle('%s/database_bnnn_%s_%s_%s.pickle' %(path, crystal_types[i], crystal_structures[i], cstypes[i]))

# Parameters for compute slcsa/atom command in LAMMPS
priors=lda_model.priors_
means=lda_model.means_    
xbar=np.dot(priors,means)

# Opération manuelle pour réduction de dimension
#data_features_reduced_test_manual[i]=np.dot(np.transpose(coefs_lda),data_features_full_test[i,:]-xbar)
  
ofile=open('dir.slcsa/mean_descriptor.dat','w')
for i in range(len(xbar)):
    ofile.write('%5.8f\n' %(xbar[i]))
ofile.close()

# On sort les coefficients de la LDA donnés par scalings_
coefs_lda = lda_model.scalings_    
ofile=open('dir.slcsa/lda_scalings.dat','w')
for i in range(len(coefs_lda)):
    ofile.write('%5.8f %5.8f %5.8f\n' %(coefs_lda[i,0], coefs_lda[i,1], coefs_lda[i,2]))
ofile.close()

# Linear Discriminant Analysis End
##################################################################

##################################################################
# Visualization purpose
fig, ax = plt.subplots(1,1,figsize=(12,9))

data=df_train_reduced.sample(5000)
label_map = dict(map(lambda i,j : (i,j) , cstypes,crystal_structures))

data['CS'] = data['CS'].map(label_map)
data=data.sort_values('CS')

sns.set_palette("rainbow", ncstypes)
a=sns.pairplot(data, diag_kind="kde", hue='CS', diag_kws={'linewidth': 1, 'alpha': 0.75}, plot_kws={'linewidth': 0, 'alpha': 0.25}, grid_kws={'despine': False},aspect=1)

sns.move_legend(a, "center",bbox_to_anchor=(0.25, 0.85), title=r'Crystal Structure', frameon=False, title_fontsize=15)

a.axes[0,0].set_yticks([])
a.axes[1,0].set_yticks([])
a.axes[2,0].set_yticks([])

a.axes[2,0].set_xticks([])
a.axes[2,1].set_xticks([])
a.axes[2,2].set_xticks([])

plt.tight_layout()
plt.savefig('LDA_separation_latent_space.png')
# Visualization purpose
##################################################################

##################################################################
# Logistic Regression on Reduced descriptor Start
print('Training a Logistic Regression Model...')

lr_model = LogisticRegression()
lr_model.fit(train_features_data_reduced, train_labels_data)

# Parameters for compute slcsa/atom command in LAMMPS
decision=lr_model.coef_
ofile=open('dir.slcsa/lr_decision.dat','w')
for i in range(len(decision)):
  for j in range(3):
    ofile.write('%5.8f ' %(decision[i,j]))
  ofile.write('\n')
ofile.close()
    
biais=lr_model.intercept_
ofile=open('dir.slcsa/lr_biais.dat','w')
for j in range(len(biais)):
  ofile.write('%5.8f ' %(biais[j]))
ofile.close()

predictions_proba_train = lr_model.predict_proba(train_features_data_reduced)
predictions_class_train = lr_model.predict(train_features_data_reduced)
print("Logistic Regression model statistics for training dataset")
print(classification_report(train_labels_data, predictions_class_train))

predictions_proba_test = lr_model.predict_proba(test_features_data_reduced)
predictions_class_test = lr_model.predict(test_features_data_reduced)
print("Logistic Regression model statistics for testing dataset")
print(classification_report(test_labels_data, predictions_class_test))

# Logistic Regression on Reduced descriptor End
##################################################################

##################################################################
# Mahalanobis distance analysis Start

def mahalanobis_distance(x,mu,icov):
  xmu=x-mu
  xmuT=np.transpose(xmu)
  left=np.dot(xmu,icov)
  right=np.dot(left,xmuT)
  return np.sqrt(right)

fig, ax = plt.subplots(1,1,figsize=(12,12))

ncomps=ncstypes-1
mean_list=[]
icov_list=[]
data_list=[]
mahafile=open('dir.slcsa/mahalanobis_file.dat','w')
for i in range(ncstypes):  
    pickle_name='bnnn_reduced_pickle_database/database_bnnn_%s_%s_%s.pickle' %(crystal_types[i], crystal_structures[i], cstypes[i])
    df_cs = pd.read_pickle(pickle_name)
    data_cs = df_cs.values[:,:ncomps]
    data_list.append(data_cs)
    meanbs=np.mean(data_cs,axis=0)
    mean_list.append(np.mean(data_cs,axis=0))
    icov=np.linalg.pinv(np.cov(np.transpose(data_cs)))
    icov_list.append(icov)
    
    # Replace with threshold
#    mahafile.write('# %s \n' %(crystal_structures[i]))
    #
#    mahafile.write('%5.4f %5.4f %5.4f\n' %(meanbs[0], meanbs[1], meanbs[2]))  
#    for m in range(3):
#        mahafile.write('%5.4f %5.4f %5.4f\n' %(icov[m,0], icov[m,1], icov[m,2]))
#    icov_list.append(np.linalg.pinv(np.cov(np.transpose(data_cs))))
#mahafile.close()

list_indexes=[[0,0],[0,1],[1,0],[1,1]]
for ind in range(ncstypes):
    data=data_list[ind]
    dmaha_to_others=np.zeros((len(data),4))
    
    for i in range(len(data)):
        xi=data[i,:]
        for j in range(ncstypes):
            muj=mean_list[j]
            icovj=icov_list[j]
            dmaha_to_others[i,j]=mahalanobis_distance(xi,muj,icovj)
    nth=1000
    err_rate_cs=np.zeros((nth,ncstypes))
    threshold_list=np.zeros(nth)
    bool_vec=np.zeros(ncstypes,dtype=bool)
    bool_vec[ind]=True
    for l in range(nth):
        threshold_list[l]=2.+0.006*l
        for k in range(ncstypes):
            if (bool_vec[k] == True):
                err_rate_cs[l,k]=100.*len(dmaha_to_others[dmaha_to_others[:,k]>threshold_list[l]])/len(data)
            else:
                err_rate_cs[l,k]=100.*len(dmaha_to_others[dmaha_to_others[:,k]<threshold_list[l]])/len(data)

    natoms_out=err_rate_cs[:,bool_vec].flatten()
    natoms_in =np.sum(err_rate_cs[:,np.invert(bool_vec)],axis=1)
    diff=np.abs(natoms_out-natoms_in)
    indmin=np.argmin(diff)    
    indi=list_indexes[ind][0]
    indj=list_indexes[ind][1]

    mahafile.write('%5.4f\n' %(threshold_list[indmin]))
    mahafile.write('%5.4f %5.4f %5.4f\n' %(mean_list[ind][0],mean_list[ind][1],mean_list[ind][2]))
    for m in range(3):
        mahafile.write('%5.4f %5.4f %5.4f\n' %(icov_list[ind][m,0], icov_list[ind][m,1], icov_list[ind][m,2]))
    
    ax.plot(threshold_list,natoms_out,label=r'$\tau_\mathrm{BCC}$',color='blue',linewidth=3.0)
    ax.plot(threshold_list,natoms_in,label=r'$\overline{\tau_\mathrm{BCC}}$',color='red',linewidth=3.0)
    ax.plot(threshold_list,diff,'--',color='green',linewidth=3.0)    
    ax.axvline(x=threshold_list[indmin],ymin=-1,ymax=101,color='black',linestyle='--',linewidth=2.0)

    ax.set_xlim([2.,8.])
    ax.set_ylim([0,10])
    ax.set_xlabel(r'Threshold d$_\mathrm{Maha}$')
    ax.set_ylabel(r'Proportion (%)')
    
    plt.tight_layout()
    plt.savefig('threshold_determination_%s_%s.png' %(crystal_types[ind], crystal_structures[ind]))
    ax.clear()
    
mahafile.close()
