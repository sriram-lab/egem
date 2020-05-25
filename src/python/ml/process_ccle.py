"""

"""

import pandas as pd

# Meta data
filename = '~/Data/Mappings/CCLE/CCLE_Summary.xlsx'
meta = pd.read_excel(filename, sheet_name='MatchedCells')
meta['CellLine'] = meta['CellLine'].str.lower()
meta['CellLine'] = meta['CellLine'].replace('-', '')

# Map CCLE_GCP Data to meta data
#filename = '~/Data/Expression/Proteomics/GCP/CCLE_GCP.csv'
#gcp = pd.read_csv(filename)
#gcp['CCL'] = gcp['CellLineName'].str.split('_').str[0]
#gcp = gcp.drop(['CellLineName', 'BroadID'], axis=1)
#gcp = gcp.set_index('CCL')

# Get match between the meta data and GCP
#merged_df = pd.merge(meta, gcp, how='outer', left_on='CellLine', right_index=True)
#merged_df.to_csv('~/Data/Expression/Proteomics/GCP/MappedGCP.csv')

#filename = '~/Data/Expression/RNASeq/CCLE/CCLE_934_TPM_EBI.tsv'
#rna = pd.read_csv(filename, delimiter='\t')
#rna.columns = rna.columns.str.split(', ').str[-1]
#rna = rna.T
#rna.columns = rna.iloc[0]
#rna = rna.drop('Genes', axis=0)
#rna.index = rna.index.str.lower()
#rna.index = rna.index.str.replace('-', '')

#merged_df = pd.merge(meta, rna, how='outer', left_on='CellLine', right_index=True)
#merged_df.to_csv('~/Data/Expression/RNASeq/CCLE/MappedRNASeq.csv')

# Map methylation data
#filename = '~/Data/Expression/Methylation/CCLE/CCLE_RRBS_tss_CpG_clusters_20181022.txt'
#methylation = pd.read_csv(filename, delimiter='\t')
#print(methylation)

# Map proteomics data
filename = '~/Data/Expression/Proteomics/CCLE/CCLE_DeepProteomics_NormalizedExpression.xlsx'
prot = pd.read_excel(filename, sheet_name='Normalized Protein Expression')
prot.columns = prot.columns.str.split('_').str[0]
print(prot.columns)