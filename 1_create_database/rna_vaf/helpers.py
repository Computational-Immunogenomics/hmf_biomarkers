import pandas as pd
import gzip

SAGE_DIR = "/data/storage/sage_vaf/"
READY_DIR = "/data/storage/ready_vaf/"

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

def filter_PASS( df: pd.DataFrame ) -> pd.DataFrame :
    return df[[True if 'PASS' in i else False for i in df['FILTER']]]

def see_annotation(vcf_path: str, match = 'INFO') -> list:
    with gzip.open(vcf_path, "rt") as ifile:
        lines = []
        for line in ifile:
            if match in line:
                lines.append(line)
        ifile.close()
    return lines

### Extract data
def dna_rna_vaf( file: str ) -> pd.DataFrame:
    a = filter_PASS(read_vcf(file))
    sample = file.split(".")[0]
    biallelic = [True if 'BIALLELIC' in i else False for i in a['INFO']]
    dna_vaf = [i.split(":")[3] for i in a[sample]]
    dna_dp = [i.split(":")[4] for i in a[sample]]
    rna_vaf = [i.split(":")[3] for i in a[sample + '_RNA\n']]
    rna_dp = [i.split(":")[4] for i in a[sample + '_RNA\n']]
    gene = [i.split("IMPACT=")[1].split(",")[0] if "IMPACT" in i else "" for i in a['INFO'].values]
    impact = [i.split("IMPACT=")[1].split(",")[2] if "IMPACT" in i else "none" for i in a['INFO']]
    purple_af = [i.split("PURPLE_AF=")[1].split(";")[0] if 'PURPLE_AF' in i else "" for i in a['INFO']]
    purple_cn = [i.split("PURPLE_CN=")[1].split(";")[0] if 'PURPLE_CN' in i else "" for i in a['INFO']]
    purple_vcn = [i.split("PURPLE_VCN=")[1].split(";")[0] if 'PURPLE_VCN' in i else "" for i in a['INFO']]
    purple_macn = [i.split("PURPLE_MACN=")[1].split(";")[0] if 'PURPLE_MACN' in i else "" for i in a['INFO']]
    subclonal = [i.split("SUBCL=")[1].split(";")[0] if 'SUBCL' in i else "" for i in a['INFO']]
    tier = [i.split("TIER=")[1].split(";")[0] if 'TIER' in i else "" for i in a['INFO']]
    mappability = [i.split("MAPPABILITY=")[1].split(";")[0] if 'MAPPABILITY' in i else "" for i in a['INFO']]
    return pd.DataFrame({'sample': sample,
                         'biallelic': biallelic,
                         'impact': impact,
                         'gene': gene,
                         'purple_af': purple_af,
                         'purple_cn': purple_cn,
                         'purple_vcn': purple_vcn,
                         'purple_macn': purple_macn,
                         'dna_vaf' : dna_vaf,
                         'rna_vaf': rna_vaf,
                         'dna_dp' : dna_dp,
                         'rna_dp' : rna_dp,
                         'subclonal': subclonal,
                         'tier' : tier,
                         'mappability' : mappability})

def filter_df( df: pd.DataFrame ) -> pd.DataFrame:
    return df[(df['gene'] != "") & (df["impact"] != "non_coding_transcript_exon_variant")]

def go(file: str) -> pd.DataFrame:
    df1 = dna_rna_vaf(file)
    return filter_df(df1)
