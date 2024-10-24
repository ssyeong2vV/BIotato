import pandas as pd
import numpy as np
import re
import sys
import warnings

vcf = sys.argv[1]
annotated = sys.argv[2]
Name = sys.argv[3]

FORMAT_HEADER = []
INFO_HEADER = []
SnpEff_HEADER = []

with open(vcf, 'r') as f :
    lines = f.readlines()
    for line in lines :
        if re.match(r"##FORMAT.+",line) :
            tmp = line.strip().split('=<ID=')
            tmp = tmp[1].split(',')
            header = tmp[0]
            FORMAT_HEADER.append(header)
        elif re.match(r"##INFO.+",line) :
            tmp = line.strip().split('=<ID=')
            tmp = tmp[1].split(',')
            header = tmp[0]
            INFO_HEADER.append(header)
            if re.match(r'##INFO=<ID=ANN,Number=.,Type=String,Description="Functional annotations:.+',line) :
                tmp = line.strip().split("'")
                header = tmp[1].split('|')
                SnpEff_HEADER.extend(header)

vcf_df = pd.read_csv(annotated, sep = '\t')

Otherinfo = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'FORMAT_DATA']
vcf_df.columns = list(vcf_df.columns)[:191] + Otherinfo

def MergeFormat(FORMAT) :
    FORMAT_merge = pd.DataFrame(columns=FORMAT_HEADER)
    for i in range(len(FORMAT)) :
        tmp = {key: '.' for key in FORMAT_HEADER}
        keys = FORMAT['FORMAT'].iloc[i].split(':')
        values = FORMAT['FORMAT_DATA'].iloc[i].split(':')
        for key, value in zip(keys, values):
            tmp[key] = value
        tmp_df = pd.DataFrame([tmp])
        FORMAT_merge = pd.concat([FORMAT_merge,tmp_df],ignore_index=True)
    return FORMAT_merge

FORMAT_df = MergeFormat(vcf_df[['FORMAT','FORMAT_DATA']])

def SplitInfo(infolist):
    tmp = {key: '.' for key in INFO_HEADER}
    for info in infolist:
        if '=' in info:
            key, value = info.split('=', 1)
            if key in tmp:
                tmp[key] = value
    return tmp
INFO_df = pd.DataFrame(columns=INFO_HEADER)
INFO_dicts = vcf_df['INFO'].apply(lambda x: SplitInfo(x.split(';')))
INFO_df = pd.DataFrame(INFO_dicts.tolist())

def SplitSnpEff(infolist) :
    tmp = {key: '.' for key in SnpEff_HEADER}
    for key,value in zip(SnpEff_HEADER,infolist) :
        tmp[key] = value
    return tmp
SnpEff_df = pd.DataFrame(columns=SnpEff_HEADER)
SnpEff_dicts = INFO_df['ANN'].apply(lambda x: SplitSnpEff(x.split('|')))
SnpEff_df = pd.DataFrame(SnpEff_dicts.tolist())

INFO_df.drop('ANN', axis=1, inplace=True)
INFO_df = pd.concat([INFO_df, SnpEff_df], axis=1)
vcf_df.drop('INFO', axis=1, inplace=True)
vcf_df.drop('FORMAT', axis=1, inplace=True)
vcf_df.drop('FORMAT_DATA', axis=1, inplace=True)
vcf_df = pd.concat([vcf_df,INFO_df,FORMAT_df], axis=1)

vcf_df['AAChange.refGeneWithVer'] = vcf_df['AAChange.refGeneWithVer'].apply(lambda x : x.split(','))
tmp_result = vcf_df.loc[np.repeat(vcf_df.index.values, vcf_df['AAChange.refGeneWithVer'].str.len())]

for i in range(len(vcf_df)) :
    if len(tmp_result['AAChange.refGeneWithVer'][i]) > 1 :
        temp = tmp_result.loc[tmp_result['AAChange.refGeneWithVer'].index.values == i,:]
        AAChange_list = temp['AAChange.refGeneWithVer'].iloc[0]
        for j in range(len(temp)) :
            temp['AAChange.refGeneWithVer'].iloc[j]=  f'{AAChange_list[j]}'
        tmp_result.loc[tmp_result['AAChange.refGeneWithVer'].index.values == i,:] = temp
tmp_result['AAChange.refGeneWithVer'] = tmp_result['AAChange.refGeneWithVer'].apply(lambda x: ''.join(x))

tmp_result[['Gene', 'NM', 'exon', 'HGVSc', 'HGVSp']] = tmp_result['AAChange.refGeneWithVer'].str.split(':', expand=True)
tmp_result.drop('AAChange.refGeneWithVer', axis=1, inplace=True)

result_header = list(tmp_result.columns)
result_header[5] = 'Func'
result_header[22] = 'gnomAD_AF'
# result_header[201] = 'Depth'
dp_columns = [i for i, col in enumerate(result_header) if col == 'DP']
result_header[dp_columns[0]] = 'Depth'

tmp_result.columns = result_header

tmp_result['AD'] = tmp_result['AD'].str.split(r'[,\|]').str[-1].astype(int)
# tmp_result['DP'] = tmp_result['DP'].replace('.', np.nan)
# tmp_result['DP'] = tmp_result['DP'].fillna(0)
tmp_result['DP'] = tmp_result['DP'].astype(int)
tmp_result['VAF.var.freq'] = tmp_result['AD']/tmp_result['DP']
tmp_result['Chr'] = 'chr'+tmp_result['Chr'].astype(str)
tmp_result['Chrom.Pos'] = tmp_result['Chr'].astype(str)+':'+tmp_result['Start'].astype(str)+'-'+tmp_result['End'].astype(str)
tmp_result.drop(labels = ['Chr','Start','End'], axis = 1,inplace = True)

tmp_result['GQ'] = tmp_result['GQ'].astype(int)
tmp_result['FS'] = tmp_result['FS'].replace('.', np.nan)
tmp_result['FS'] = tmp_result['FS'].astype(float)
tmp_result['gnomAD_AF'] = tmp_result['gnomAD_AF'].replace('.', np.nan)
tmp_result['gnomAD_AF'] = tmp_result['gnomAD_AF'].astype(float)
tmp_result['select'] = ""


tmp_result.columns = tmp_result.columns.str.strip()
tmp_result['Rank'] = tmp_result['Rank'].str.split('/').str[-1]
tmp_result['exon'] = tmp_result['exon'].str.extract('(\d+)')
tmp_result['exon'] = tmp_result['exon'].astype(str) + '_' + tmp_result['Rank']
tmp_result.loc[tmp_result['Func'] == 'intronic', 'exon'] = '.'
tmp_result.drop(columns=['Rank'], inplace=True)

if 'FS' not in tmp_result.columns:
    tmp_result = tmp_result.assign(FS='')

selected = []
with open('/labmed/01.ALL/03.python/lsy/00.CODE/hg19_RefSeq_Selected.txt', 'r') as file :
    lines = file.readlines()
    for line in lines :
        l = line.strip()
        selected.append(l)
tmp_result['main'] = tmp_result['NM'].apply(lambda x : '     O' if x in selected else "")

front_columns = ['select','main','CLNSIG', 'Chrom.Pos', 'Gene', 'NM', 'HGVSc', 'HGVSp', 'exon', 'Func', 
    'Annotation', 'Annotation_Impact', 'VAF.var.freq', 'AD', 'DP', 'GQ', 'FS', 'FILTER', 
    'eQTLGen_snp_id', 'LOF', 'NMD', 'ONC', 'ONCDN', 'ONCDISDB', 'ONCREVSTAT', 'CLNALLELEID', 
    'CLNDN', 'CLNDISDB', 'CLNREVSTAT']
back_columns = ['Ref', 'Alt','Gene.refGeneWithVer', 'GeneDetail.refGeneWithVer', 'ExonicFunc.refGeneWithVer']
middle_columns = [col for col in tmp_result.columns if col not in front_columns + back_columns]
new_column_order = front_columns + middle_columns + back_columns
result = tmp_result[new_column_order]

result.to_excel(f'{Name}.xlsx', index = False)