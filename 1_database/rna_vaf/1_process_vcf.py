#!/usr/bin/env python
# coding: utf-8

# In[10]:


import os
import pandas as pd
import gzip
from helpers import *


# ### Newer helpers

# - Newer

# In[11]:


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
            df = metadata.join(info).join(vaf)
            df = df[(df['gene'].notna()) & (df['FILTER'] == "PASS")] 
            df['sampleId'] = sample
            df = df.drop(['FILTER', 'impact'], axis=1)
        return df
    except Exception as e:
        print(f"Unexpected error: {e}")
        return None                 


# #### 0 - Go!

# In[3]:


os.chdir(SAGE_DIR)


# In[4]:


vafs = [i for i in os.listdir() if 'purple.somatic' in i]


# In[5]:


vaf_dfs = []; i = 0
for f in vafs:
    print("Processing " + f)
    i = i+1
    print(i)
    vaf_dfs.append( vamos(f) )


# #### 1 - Pre-filtering of the data

# In[6]:


tmp = pd.concat(vaf_dfs)


# #### 2 - Send it!

# In[7]:


tmp.to_csv( READY_DIR + "vaf_ready_test.csv")

