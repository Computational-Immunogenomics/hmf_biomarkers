#!/usr/bin/env python
# coding: utf-8

# In[1]:


import sys
import os
import gzip
import pandas as pd 
import warnings
warnings.filterwarnings('ignore')


# In[2]:


sys.path.insert(0, f'{os.path.dirname(os.getcwd())}')


# In[3]:


from map import I_DIR, SOMATIC_DIR, O_DIR
#SOMATIC_DIR = I_DIR + 'somatics/'
#O_DIR = I_DIR + 'biomarkers/database/'


# ### Helper Functions 

# In[5]:


### Functions to Read VCF file ###
def fps(sample: str) -> str:
    return SOMATIC_DIR + sample + "/purple/" + sample + ".purple.somatic.vcf.gz"

def read_vcf(vcf_path: str) -> pd.DataFrame:
    names = get_vcf_col_names(vcf_path)
    return pd.read_csv(vcf_path, compression='gzip', comment='#', chunksize=10000, delim_whitespace=True, header=None, names=names).read()

def get_vcf_col_names(vcf_path: str) -> list:
    with gzip.open(vcf_path, "rt") as ifile:
          for line in ifile:
            if line.startswith("#CHROM"):
                  vcf_names = [x for x in line.split('\t')]
                  break
    ifile.close()
    return vcf_names

### Functions to extract essential fields for feature creation ###
def get_metadata_df(vcf: pd.DataFrame) -> pd.DataFrame:
    return vcf[["#CHROM", "POS",  "REF", "ALT", "FILTER"]].rename(columns = {'#CHROM':'chromosome', "POS" : "position"})

def filter_PASS( df: pd.DataFrame ) -> pd.DataFrame :
    return df[[True if 'PASS' in i else False for i in df['FILTER']]]

#### pass INFO column #### 
def get_impact_fields( s: str ) :
    impact_fields = s.split("IMPACT=")[1].split(";")[0].split(",")
    return impact_fields    
        
def get_Gene(s: str):
    if 'IMPACT=' in s: 
        return get_impact_fields( s )[0]
    else:
        return pd.NA 

def get_Transcript(s: str):
    if 'IMPACT=' in s: 
        return get_impact_fields( s )[1]
    else:
        return pd.NA 

def get_Label(s: str):
    if 'IMPACT=' in s: 
        return get_impact_fields( s )[2]
    else:
        return pd.NA    

def get_field(s: str, f: str, n: bool) -> bool:
    if f in s : 
        field = s.split(f)[-1].split(';')[0]
        if n:
            return float(field)
        else: 
            return field
    else:
        return pd.NA

def get_biallelic(s: str) -> bool:
    if 'BIALLELIC' in s : 
        return True
    else:
        return False

def extract_info_fields(s: str) -> dict:
    return {
        'gene': get_Gene(s), 
        'transcript' : get_Transcript(s), 
        'annotation': get_Label(s),
        'purple_af': get_field(s, "PURPLE_AF=", True),
        'purple_cn': get_field(s, "PURPLE_CN=", True),
        'purple_vcn': get_field(s, "PURPLE_VCN=", True),
        'purple_macn': get_field(s, "PURPLE_MACN=", True),
        'tier': get_field(s, "TIER=", False),
        'tnc' : get_field(s, "TNC=", False),        
        'biallelic': get_biallelic(s),
        'subclonal': get_field(s, "SUBCL=", True),
     }    

def get_info_df(vcf: pd.Series) -> pd.DataFrame:
    return pd.DataFrame([extract_info_fields(i) for i in vcf['INFO']])
    
def vamos(sampleId: str) -> pd.DataFrame:
    try:
        fp = fps(sampleId)
        if os.path.exists(fp):
            vcf = read_vcf(fp)
            metadata = get_metadata_df( vcf )
            info = get_info_df( vcf )
            ready = get_metadata_df(vcf).join(get_info_df(vcf))
            ready = ready.query('FILTER=="PASS" & annotation == annotation & gene != ""')
            #ready = ready.query('FILTER=="PASS" & (annotation == annotation | tier == "HOTSPOT") & (annotation not in ("3_prime_UTR_variant","5_prime_UTR_variant","intron_variant","upstream_gene_variant","non_coding_transcript_exon_variant"))')    
            ready[['sampleId']] = sampleId
        return ready
    except Exception as e:
        print(f"Unexpected error: {e}")
        return None                 


# ### Get file paths 

# In[6]:


files = os.listdir(SOMATIC_DIR)
run_files = [i for i in files if "." not in i]


# ### Run and save output

# In[7]:


chunk_size = 702


# In[9]:


dfs = []; k = 0; j = 0
for i in run_files:
    k = k + 1
    print(k); print(i)
    dfs.append(vamos(i))
    if (k % chunk_size) == 0:
        j = j + 1
        pd.concat(dfs).to_csv(O_DIR + 'somatic_exome/somatic_exome_chunk_' + str(j) + '.csv', index = False)
        dfs = []
    if k == len(run_files):
        j = j + 1
        pd.concat(dfs).to_csv(O_DIR + 'somatic_exome/somatic_exome_chunk_' + str(j) + '.csv', index = False)    

