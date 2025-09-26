#!/usr/bin/env python
# coding: utf-8

# In[55]:


import os
import pandas as pd
import gzip
from helpers import *


# ### Helpers 

# - Old

# - Newer

# In[43]:


### Functions to Read VCF file ###
def get_vcf_col_names(vcf_path: str) -> list:
    with gzip.open(vcf_path, "rt") as ifile:
          for line in ifile:
            if line.startswith("#CHROM"):
                  vcf_names = [x for x in line.split('\t')]
                  break
    ifile.close()
    return vcf_names

def read_vcf(vcf_path: str) -> pd.DataFrame:
    names = get_vcf_col_names(vcf_path)
    return pd.read_csv(vcf_path, compression='gzip', comment='#', chunksize=10000, sep='\s+', header=None, names=names).read()

### Metadata ###
def get_metadata_df(vcf: pd.DataFrame) -> pd.DataFrame:
    return vcf[["#CHROM", "POS",  "REF", "ALT", "FILTER"]].rename(columns = {'#CHROM':'chromosome', "POS" : "position"})

#### pass INFO column #### 
def get_impact_fields( s: str ) :
    impact_fields = s.split("IMPACT=")[1].split(";")[0].split(",")
    return impact_fields    
        
def get_impact_field(s: str, field: str):
    if 'IMPACT=' in s: 
        impact = get_impact_fields( s )
        if field == 'gene':
            return impact[0]
        elif field == 'transcript':
            return impact[1]
        elif field == 'impact':
            return impact[2]
    else:
        return pd.NA 
        
def extract_info_fields(s: str) -> dict:
    return {
        'gene': get_impact_field(s, "gene"), 
        'transcript' : get_impact_field(s, "transcript"), 
        'impact': get_impact_field(s, "impact")
     }    

def get_info_df(vcf: pd.Series) -> pd.DataFrame:
    return pd.DataFrame([extract_info_fields(i) for i in vcf['INFO']])

### Extra RNA and DNA depth/vaf ### 
def get_rna_dna_dp_vaf(sample: str, vcf: pd.DataFrame) -> pd.DataFrame:
    return pd.DataFrame(
    {'dna_vaf' : [i.split(":")[3] for i in vcf[sample]],
     'dna_dp' : [i.split(":")[4] for i in vcf[sample]],
     'rna_vaf' : [i.split(":")[3] for i in vcf[sample + '_RNA\n']],
     'rna_dp' : [i.split(":")[4] for i in vcf[sample + '_RNA\n']]
    })

### Together ### 
def vamos(fp: str) -> pd.DataFrame:
    sample = fp.split(".")[0]
    try:
        if os.path.exists(fp):
            vcf = read_vcf(fp)
            metadata = get_metadata_df( vcf )
            info = get_info_df( vcf )
            vaf = get_rna_dna_dp_vaf( sample, vcf)
            together = metadata.join(info).join(vaf)
            ready = together.query('FILTER=="PASS" & gene != "")')    
            ready[['sampleId']] = sample
        return ready
    except Exception as e:
        print(f"Unexpected error: {e}")
        return None                 


# In[46]:


fp = vafs[0]


# In[54]:


together


# In[65]:


#vamos(fp)
sample = fp.split(".")[0]
vcf = read_vcf(fp)
metadata = get_metadata_df( vcf )
info = get_info_df( vcf )
vaf = get_rna_dna_dp_vaf( sample, vcf)
df = metadata.join(info).join(vaf)

df_go = df[(df['gene'].notna()) & (df['FILTER'] == "PASS")]


# In[67]:


df_go[df_go['impact'] == "missense_variant"]


# #### 0 - Go!

# In[10]:


os.chdir(SAGE_DIR)


# In[45]:


vafs = [i for i in os.listdir() if 'purple.somatic' in i]


# In[297]:


vaf_dfs = []
for f in vafs[:3]:
    print("Processing " + f)
    vaf_dfs.append( go(f) )


# #### 1 - Pre-filtering of the data

# In[298]:


tmp = pd.concat(vaf_dfs)


# In[299]:


tmp['rna_dp'] = tmp['rna_dp'].astype(int)


# In[300]:


tmp.to_csv(READY_DIR + "vaf_test.csv")


# In[301]:


READY_DIR + "vaf_test.csv"


# In[52]:


def filter_by_column(df, column, value):
    return df[df[column] > value]


# In[53]:


filtered_df = tmp.pipe(filter_by_column, 'rna_dp', 4)


# #### 2 - Send it!

# In[54]:


filtered_df.to_csv( READY_DIR + "vaf_ready.csv")

