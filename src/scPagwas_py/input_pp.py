import scanpy as sc
import numpy as np
import pandas as pd
from tqdm import tqdm
from scipy.sparse import csc_matrix
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from sklearn.utils.extmath import randomized_svd
import pybedtools
from pybedtools import BedTool

def scdata_process(Pagwas=None,adata=None, cell_type_col='annotation'):
    # 根据细胞类型计算平均表达量
    cell_types = adata.obs[cell_type_col].unique()
    avg_expr_matrix = {}
    for cell_type in cell_types:
        cells_of_type = adata[adata.obs[cell_type_col] == cell_type]
        avg_expr = cells_of_type.X.mean(axis=0)
        avg_expr_matrix[cell_type] = avg_expr.A1 if isinstance(avg_expr, np.matrix) else avg_expr
    # 将平均表达量矩阵转换为DataFrame
    avg_expr_matrix = pd.DataFrame(avg_expr_matrix, index=adata.var_names)
    if Pagwas is None:
        Pagwas = {}
    Pagwas['avg_expr_matrix'] = avg_expr_matrix
    adata.obs = adata.obs.rename(columns={cell_type_col: 'annotation'})
    Pagwas['Celltype_anno'] = adata.obs
    X = adata.X
    # 2. 获取基因名和细胞名
    gene_names = adata.var_names  # 获取基因名称
    cell_names = adata.obs_names  # 获取细胞名称
    if hasattr(X, 'toarray'):
        X_dense = X.toarray()
    else:
        X_dense = X  # 如果已经是密集矩阵则直接使用
    X_dense=X_dense.T
    # 创建 DataFrame
    df = pd.DataFrame(X_dense, index=gene_names, columns=cell_names)
    Pagwas['data_mat'] = df
    Pagwas['Celltype_anno']['cellnames'] = Pagwas['Celltype_anno'].index
    return Pagwas


    """
    读取并预处理GWAS数据。
    参数:
    Pagwas : dict
        包含GWAS数据的字典
    gwas_data : pandas.DataFrame
        GWAS数据框
    maf_filter : float
        次要等位基因频率（MAF）的过滤阈值
    Sex_filter : bool
        是否过滤性染色体
    MHC_filter : bool
        是否过滤主要组织相容性复合体（MHC）区域
    gwas_z_filter : float
        |z|值的过滤阈值
    返回:
    Pagwas : dict
        更新后的包含预处理过的GWAS数据的字典
    """
from pybedtools import BedTool

def snp_to_gene(gwas_data, block_annotation, marg=10000):
    gwas_data['pos_start'] = gwas_data['pos']
    gwas_data['pos_end'] = gwas_data['pos']
    gwas_bed = BedTool.from_dataframe(gwas_data[['chrom', 'pos_start', 'pos_end', 'rsid']]).sort()
    block_annotation['start'] = block_annotation['start'] + 1 - marg
    block_annotation['start'] = block_annotation['start'].apply(lambda x: max(x, 0))
    block_annotation['end'] = block_annotation['end'] + marg
    block_bed = BedTool.from_dataframe(block_annotation[['chrom', 'start', 'end', 'label']]).sort()
    nearest = gwas_bed.closest(block_bed, d=True)
    nearest_df = nearest.to_dataframe(names=['chrom', 'pos_start', 'pos_end', 'rsid', 
                                             'chrom_gene', 'gene_start', 'gene_end', 
                                             'label', 'distance'])
    nearest_df = nearest_df[nearest_df['distance'] == 0]
    map_result = nearest_df[['rsid', 'label', 'pos_start', 'distance']]
    map_result.columns = ['rsid', 'label', 'pos', 'Distance']
    # 保留rsid列，并同时将其设为索引
    # 使用 .loc 进行赋值，避免 SettingWithCopyWarning
    map_result.loc[:, 'rsid_index'] = map_result['rsid']
    map_result.set_index('rsid_index', inplace=True)

    return map_result


def GWAS_summary_input(Pagwas=None, gwas_data=None, maf_filter=0.01, 
                       Sex_filter=True, MHC_filter=True, gwas_z_filter=-1,
                       block_annotation=None, marg=10000):
    if gwas_data is None or not isinstance(gwas_data, pd.DataFrame):
        raise ValueError("The input data is not in dataframe format!")
    necessary_cols = ["chrom", "pos", "rsid", "beta", "se", "maf"]
    if not all(col in gwas_data.columns for col in necessary_cols):
        raise ValueError("Necessary columns are missing: chrom, pos, beta, se, maf")
    gwas_data['pos'] = pd.to_numeric(gwas_data['pos'], errors='coerce')
    gwas_z_filter = float(gwas_z_filter)
    maf_filter = float(maf_filter)
    gwas_data['chrom'] = gwas_data['chrom'].astype(str)
    if not gwas_data['chrom'].str.startswith('chr').any():
        gwas_data['chrom'] = 'chr' + gwas_data['chrom']
    if 'maf' in gwas_data.columns:
        gwas_data['maf'] = np.where(gwas_data['maf'] > 0.5, 1 - gwas_data['maf'], gwas_data['maf'])
        gwas_data = gwas_data[gwas_data['maf'] > maf_filter]
    gwas_data = gwas_data.drop_duplicates(subset=['rsid'])
    if Sex_filter:
        gwas_data = gwas_data[~gwas_data['chrom'].isin(['chr23', 'chrX', 'chrY'])]
    if MHC_filter:
        mhc_region = gwas_data[(gwas_data['chrom'] == 'chr6') & (gwas_data['pos'] > 25000000) & (gwas_data['pos'] < 34000000)]
        gwas_data = gwas_data[gwas_data['chrom'] != 'chr6']
        gwas_data = pd.concat([gwas_data, mhc_region])
    if gwas_z_filter > 0:
        gwas_data['abs_z'] = (gwas_data['beta'] / gwas_data['se']).abs()
        gwas_data = gwas_data[gwas_data['abs_z'] < gwas_z_filter]
    gwas_data = gwas_data[gwas_data['beta'] < 1]
    if Pagwas is None:
        Pagwas = {}
    Pagwas['gwas_data'] = gwas_data
    Pagwas['snp_gene_df'] = snp_to_gene(gwas_data=gwas_data, block_annotation=block_annotation, marg=marg)
    return Pagwas

# 示例使用
# 假设gwas_data是一个pandas DataFrame，包含必要的列
# gwas_data = pd.read_csv('path_to_your_gwas_data.csv')
# Pagwas = GWAS_summary_input(Pagwas={}, gwas_data=gwas_data)
# print(Pagwas['gwas_data'])



import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from scipy.sparse import csr_matrix

def pathway_pcascore_run(Pagwas=None, Pathway_list=None, min_pathway_size=10, max_pathway_size=1000):
    if Pathway_list is None:
        raise ValueError("Pathway_list data not loaded")
    if not isinstance(Pathway_list, dict):
        raise ValueError("Pathway_list data must be a dictionary")
    # Filter pathways by size
    Pathway_list = {k: v for k, v in Pathway_list.items() if min_pathway_size <= len(v) <= max_pathway_size}
    
    # Validate gene overlap
    pa_gene = set().union(*Pathway_list.values())
    if len(set(Pagwas['data_mat'].index) & pa_gene) < len(pa_gene) * 0.1:
        raise ValueError("Less than 10% intersecting genes between Single_data and Pathway_list. Check gene names.")
    # Filter out scarce pathways
    Pathway_list = {k: v for k, v in Pathway_list.items() if len(set(v) & set(Pagwas['data_mat'].index)) > 2}
    Pagwas['rawPathway_list'] = Pathway_list
    # Filter out genes with no expression in single cells for each pathway
    celltypes = Pagwas['Celltype_anno']['annotation'].unique()
    pana_list = []
    for celltype in celltypes:
        scCounts = Pagwas['data_mat'].loc[:, Pagwas['Celltype_anno']['cellnames'][Pagwas['Celltype_anno']['annotation'] == celltype]]
        scCounts = scCounts.loc[scCounts.sum(axis=1) != 0, :]
        proper_gene_names = scCounts.index
        pana = [k for k, v in Pathway_list.items() if len(set(v) & set(proper_gene_names)) > 2]
        pana_list.append(pana)
    valid_pathways = set.intersection(*map(set, pana_list))
    Pathway_list = {k: v for k, v in Pathway_list.items() if k in valid_pathways}
    Pagwas['Pathway_list'] = Pathway_list
    Pagwas['rawPathway_list'] = Pathway_list.copy()
    print("* Start to get Pathway SVD score!")
    scPCAscore_list = []
    for celltype in celltypes:
        scCounts = Pagwas['data_mat'].loc[:, Pagwas['Celltype_anno']['cellnames'][Pagwas['Celltype_anno']['annotation'] == celltype]]
        scPCAscore = pathway_pca_test(Pathway_list=Pathway_list, scCounts=scCounts)
        print(celltype)
        scPCAscore_list.append(scPCAscore)
    # Remove pathways with issues
    pa_remove = set().union(*[x[2] for x in scPCAscore_list])
    # Merge PCA scores
    pca_df = []
    for i, (df, _, _) in enumerate(scPCAscore_list):
        df = df[~df['name'].isin(pa_remove)]
        df['celltype'] = celltypes[i]
        pca_df.append(df)
    pca_df = pd.concat(pca_df)
    pca_scoremat = pca_df.pivot(index='name', columns='celltype', values='score')
    pca_scCell_mat = pd.concat([x[1].loc[~x[1].index.isin(pa_remove)] for x in scPCAscore_list], axis=1)
    Pagwas['pca_scCell_mat'] = DataMatrix(pca_scCell_mat)
    Pagwas['pca_cell_df'] = DataMatrix(pca_scoremat)
    if pa_remove:
        Pagwas['Pathway_list'] = {k: v for k, v in Pagwas['Pathway_list'].items() if k not in pa_remove}
        Pagwas['rawPathway_list'] = {k: v for k, v in Pagwas['rawPathway_list'].items() if k not in pa_remove}
    
    return Pagwas

def DataMatrix(df):
    """
    Converts a pandas DataFrame into a dictionary containing:
    - The numpy array of the data
    - The index (row labels)
    - The columns (column labels)
    
    :param df: pandas DataFrame
    :return: Dictionary containing 'data', 'index', and 'columns'
    """
    return {
        'data': df.to_numpy(),  # Convert to numpy array
        'index': df.index,      # Row labels (index)
        'columns': df.columns   # Column labels
    }
        
import pandas as pd
from sklearn.decomposition import PCA
import numpy as np

def pathway_pca_test(Pathway_list, scCounts):
    if scCounts.index.duplicated().any():
        raise ValueError("Duplicate gene names are not allowed - please reduce")
    if scCounts.columns.duplicated().any():
        raise ValueError("Duplicate cell names are not allowed - please reduce")
    if scCounts.index.isna().any():
        raise ValueError("NA gene names are not allowed - please fix")
    if scCounts.columns.isna().any():
        raise ValueError("NA cell names are not allowed - please fix")
    scCounts = scCounts.T
    cm = scCounts.mean(axis=0)
    proper_gene_names = scCounts.columns
    papca = {}
    pa_remove = []
    for k, Pa_id in Pathway_list.items():
        lab = proper_gene_names.isin(set(Pa_id) & set(proper_gene_names))
        mat = scCounts.loc[:, lab]
        try:
            pca = PCA(n_components=1)
            pca.fit(mat)
            pcs = pca.transform(mat)
            # Adjust scores with correlation-based sign
            pcs_scores = pcs[:, 0] * np.sign(np.corrcoef(pcs[:, 0], mat.mean(axis=1))[0, 1])
            papca[k] = {'scores': pcs_scores, 'rotation': pca.components_[0]}
        except Exception as e:
            pa_remove.append(k)
    # Construct the DataFrame for variance scores
    vdf = pd.DataFrame({
        'name': list(papca.keys()),
        'score': [np.var(p['scores']) for p in papca.values()]
    })
    # Construct the DataFrame for the PCA scores
    vscore = pd.DataFrame({
        k: p['scores'] for k, p in papca.items()
    }).T
    return vdf, vscore, pa_remove


def pathway_annotation_input(Pagwas, block_annotation):
    # Check the input types
    if not isinstance(Pagwas['Pathway_list'], dict):
        raise ValueError("Pathway_list must be a dictionary")

    if not isinstance(block_annotation, pd.DataFrame):
        raise ValueError("block_annotation must be a pandas DataFrame")

    required_columns = {'chrom', 'start', 'end', 'label'}
    if not required_columns.issubset(block_annotation.columns):
        raise ValueError("block_annotation must contain columns: 'chrom', 'start', 'end', 'label'")

    # Remove duplicated labels in block_annotation
    block_annotation = block_annotation.drop_duplicates(subset='label')

    # Intersect variable features with gene labels
    proper_gene_names = set(Pagwas['data_mat'].index).intersection(Pagwas['snp_gene_df']['label'])

    if len(set().union(*Pagwas['Pathway_list'].values()).intersection(proper_gene_names)) < 1:
        raise ValueError("No match between Pathway genes and VariableFeatures")

    # Filter the Pathway_list to only include valid genes
    Pathway_list = {
        pa: list(set(genes).intersection(proper_gene_names))
        for pa, genes in Pagwas['Pathway_list'].items()
    }

    # Filter out empty pathways
    Pathway_list = {k: v for k, v in Pathway_list.items() if len(v) > 0}

    # Keep only pathways present in pca_cell_df's rows
    pca_cell_df = pd.DataFrame(Pagwas['pca_cell_df']['data'], index=Pagwas['Pathway_list'].keys())

    # 然后再执行集合交集操作
    Pa_index = set(Pathway_list.keys()).intersection(pca_cell_df.index)
    Pathway_list = {k: Pathway_list[k] for k in Pa_index}

    print("Obtaining pathway block information...")

    # Filter pathways based on block_annotation
    valid_pa_index = [
        pa for pa in Pathway_list if block_annotation['label'].isin(Pathway_list[pa]).sum() >= 10
    ]
    Pathway_list = {pa: Pathway_list[pa] for pa in valid_pa_index}

    # Create a DataFrame for each pathway and label it
    paan_df = []
    for pa, genes in Pathway_list.items():
        pa_block = block_annotation[block_annotation['label'].isin(genes)].copy()
        pa_block['pathway'] = pa
        paan_df.append(pa_block)

    # Sort the blocks by chromosome and start position
    pathway_blocks = {
        pa: df.sort_values(['chrom', 'start']) for pa, df in zip(Pathway_list.keys(), paan_df)
    }

    # Update Pagwas object with filtered Pathway_list and pathway_blocks
    Pagwas['Pathway_list'] = Pathway_list
    Pagwas['pathway_blocks'] = pathway_blocks

    return Pagwas



