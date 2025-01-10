#!/home/mmoro/SPATIAL/Cell_neigh/Neigh_v2/SCOTIA/scotia/bin/python
import numpy as np
import pandas as pd
import seaborn as sns
import scotia
import matplotlib.pyplot as plt
from matplotlib import cm
from scipy.spatial import distance_matrix
from scipy.stats import ranksums
from matplotlib.collections import LineCollection
from statsmodels.stats.multitest import multipletests
import warnings
warnings.filterwarnings("ignore")
import sys

def read_data(file):
    try:
        df = pd.read_csv(file, header = 0, index_col=None, sep = "\t")
    except:
        df = pd.DataFrame()
    return df
  

known_lr_pairs = pd.read_csv('/home/mmoro/SPATIAL/Mackensy_analysis/CosMx_Prot/3_SCOTIA/Files/lr_pair.csv',sep=',', index_col = None)
print(known_lr_pairs.head())

samples = [f"patient{i}" for i in range(1, 34)]

for patient in samples:
    
  meta_df = pd.read_csv('/home/mmoro/SPATIAL/Mackensy_analysis/CosMx_Prot/3_SCOTIA/Files/'+patient+'_meta_def.csv',sep=',',index_col=None)
  print(meta_df.head())
  exp_df = pd.read_csv('/home/mmoro/SPATIAL/Mackensy_analysis/CosMx_Prot/3_SCOTIA/Files/'+patient+'_exp_def.csv',sep=',',index_col=None)
  print(exp_df.head())
  
  #Obtain all FOVS
  fovs = meta_df.iloc[:,1].unique()
  
  #get DBSCAN cell clusters
  #for fov in set(meta_df['fov']):
  for fov in fovs:
      print(fov)
      celltype = []
      cell_idx = []
      meta_df_fov = meta_df[meta_df['fov']==fov]
      meta_df_fov['index'] = range(meta_df_fov.shape[0])
      meta_df_fov.to_csv('/home/mmoro/SPATIAL/Mackensy_analysis/CosMx_Prot/3_SCOTIA/Files/'+patient+'_fov_'+str(fov)+".csv",header = True, index = False, sep = "\t")
  
      cell_type_l = list(set(meta_df_fov['annotation']))
      
      for ct in cell_type_l:
          ###cluster cells
          #print(ct)
          meta_df_sel = meta_df_fov[meta_df_fov['annotation']==ct]
          X = np.array(meta_df_sel[['x_positions','y_positions']])
          if X.shape[0] >= 5:
              print(ct)
              idx_l, fi_eps = scotia.dbscan_ff_cell(X, X_index_arr=np.array(meta_df_sel['index']),min_cluster_size=5,eps_range = list(range(15,300,5)))
              if len(idx_l)>0:
                  celltype += [ct for x in idx_l]
                  cell_idx += idx_l
                  
                  
                  
      tmp_df = pd.DataFrame([celltype,cell_idx]).T
      np.save('/home/mmoro/SPATIAL/Mackensy_analysis/CosMx_Prot/3_SCOTIA/Files/'+patient+'_fov_'+str(fov)+'_dbscan.cell.clusters.test',tmp_df)
        
  #ot 
  #gene expression normalization factor
  exp_df_all = exp_df
  exp_df_norm = exp_df_all.iloc[:,2:]
  exp_df_norm = exp_df_norm[exp_df_norm>0]
  df_quantile = exp_df_norm.quantile(q=0.99,axis = 0) 
  
  for fov in set(meta_df['fov']):
  #for fov in [9]:
      print(fov)
      #clustering result
      cluster_df = pd.DataFrame(np.load('/home/mmoro/SPATIAL/Mackensy_analysis/CosMx_Prot/3_SCOTIA/Files/'+patient+'_fov_'+str(fov)+'_dbscan.cell.clusters.test.npy',allow_pickle=True))
      cluster_df.columns = ['cell_type','cell_idx']
      
      #coordinates
      meta_df_fov = meta_df[meta_df['fov']==fov]
      meta_df_fov.index = range(meta_df_fov.shape[0])
      cell_id_all = np.array(range(meta_df_fov.shape[0]))
      coord = np.array(meta_df_fov[['x_positions','y_positions']])
      S_all_arr = distance_matrix(coord,coord)
      
      #expression
      exp_df_fov = exp_df[exp_df['fov']==fov].iloc[:,2:]
      exp_df_fov = exp_df_fov/df_quantile
      exp_df_fov[exp_df_fov>1]=1
      exp_df_fov.index = cell_id_all
  
      #select potentially communicating cell cluster pairs (spatially adjacent)
      S_all_arr_new = scotia.sel_pot_inter_cluster_pairs(S_all_arr,cluster_df)
  
      #optimal transport between source and target cells
      ga_df_final = scotia.source_target_ot(S_all_arr_new, exp_df_fov, meta_df_fov, known_lr_pairs)
      
      
      if ga_df_final.shape[0]>0:
          ga_df_final.columns = ['source_cell_idx','receptor_cell_idx','likelihood','ligand_recptor','source_cell_type','target_cell_type']
          ga_df_final.to_csv('/home/mmoro/SPATIAL/Mackensy_analysis/CosMx_Prot/3_SCOTIA/Files/'+patient+'_fov_'+str(fov)+".ot.csv",header = True, index = False, sep = "\t")
  
          #post-processing of ot results by calculating averaged likelihoods
          ga_df_final['cell_pairs'] = ga_df_final['source_cell_type']+"_"+ga_df_final['target_cell_type']
          ga_df_final.to_csv('/home/mmoro/SPATIAL/Mackensy_analysis/CosMx_Prot/3_SCOTIA/Files/'+patient+'_fov_'+str(fov)+".ot.csv",header = True, index = False, sep = "\t")
          # final_summary = scotia.post_ot(ga_df_final,label = patient)
          # final_summary.to_csv('/home/mmoro/SPATIAL/Mackensy_analysis/CosMx_Prot/3_SCOTIA/Files/'+patient+'_fov_'+str(fov)+".post.ot.csv",header = True, index = False, sep = "\t")
      
