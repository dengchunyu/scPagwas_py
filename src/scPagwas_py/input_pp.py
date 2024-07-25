import scanpy as sc
import numpy as np
import pandas as pd
from tqdm import tqdm
from scipy.sparse import csc_matrix
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from sklearn.utils.extmath import randomized_svd

def scdata_process(Pagwas=None,adata=None, cell_type_col='cell_type'):
    """
    使用scanpy预处理单细胞数据并计算每个细胞类型的平均表达量。
    
    参数:
    Pagwas : dict
    包含GWAS数据的字典
    adata : AnnData
        原始单细胞表达数据
    cell_type_col : str
        细胞类型注释信息列的名称
    is_normalized : bool
        数据是否已经标准化处理
    is_filtered : bool
        数据是否已经经过过滤低质量细胞
    
    返回:
    avg_expr_matrix : pandas.DataFrame
        每个细胞类型的平均表达矩阵
    """
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
    cell_names = adata.obs['cellnames']  # 获取细胞名称
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
def GWAS_summary_input(Pagwas=None, gwas_data=None, maf_filter=0.01, 
                       Sex_filter=True, MHC_filter=True, gwas_z_filter=-1):
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
    return Pagwas

# 示例使用
# 假设gwas_data是一个pandas DataFrame，包含必要的列
# gwas_data = pd.read_csv('path_to_your_gwas_data.csv')
# Pagwas = GWAS_summary_input(Pagwas={}, gwas_data=gwas_data)
# print(Pagwas['gwas_data'])


def pathway_pca_test(pathway_list, sc_counts):
    """
    计算每个路径的PCA分数。

    参数:
    pathway_list : dict
        基因列表，键是路径名称，值是基因符号列表。
    sc_counts : pandas.DataFrame
        单细胞计数矩阵。

    返回:
    list
        包含PCA分数的数据帧和PCA分数矩阵。
    """
    # 转置sc_counts
    sc_counts = sc_counts.T
    cm = np.array(sc_counts.mean(axis=0))
    proper_gene_names = sc_counts.columns

    # 计算每个路径的PCA
    papca = []
    for pathway_name, genes in pathway_list.items():
        lab = proper_gene_names.isin(genes)
        mat = sc_counts.iloc[:, lab]
        try:
            u, s, vt = randomized_svd(mat, n_components=1, random_state=0)
            scores = u @ np.diag(s)
            scores = scores - cm[lab] @ vt.T
            
            cs = np.array([np.sign(np.corrcoef(scores[:, i], (mat * np.abs(vt[i, :])).mean(axis=1))[0, 1]) for i in range(scores.shape[1])])
            scores *= cs
            vt *= cs[:, np.newaxis]
            
            papca.append({
                'xp': {
                    'd': s / np.sqrt(mat.shape[0]),
                    'rotation': vt.T,
                    'scores': scores
                },
                'n': mat.shape[1]
            })
        except Exception as e:
            papca.append(None)
    pa_remove = [name for name, result in zip(pathway_list.keys(), papca) if result is None]
    papca = [result for result in papca if result is not None]
    vdf = pd.DataFrame(
        np.concatenate(
            [np.column_stack((np.repeat(i, len(pca['xp']['d'])), pca['xp']['d'], np.repeat(pca['n'], len(pca['xp']['d'])), np.arange(1, len(pca['xp']['d']) + 1))) for i, pca in enumerate(papca)]
        ),
        columns=['i', 'var', 'n', 'npc']
    )
    
    vscore = pd.DataFrame(
        np.concatenate([pca['xp']['scores'] for pca in papca]),
        index=[name for name, result in zip(pathway_list.keys(), papca) if result is not None],
        columns=sc_counts.index
    )
    
    n_cells = sc_counts.shape[0]
    vdf['exp'] = 1  # 简化处理，R代码中使用了RMTstat::qWishartMax函数，这里用1代替
    vdf['var'] = vdf['var'] / vdf['exp']
    
    df = pd.DataFrame({
        'name': [name for name, result in zip(pathway_list.keys(), papca) if result is not None],
        'score': vdf['var']
    })
    
    return [df, vscore, pa_remove]

# 示例使用
# pathway_list = {'pathway1': ['gene1', 'gene2'], 'pathway2': ['gene3', 'gene4']}
# sc_counts = pd.DataFrame(np.random.rand(100, 4), columns=['gene1', 'gene2', 'gene3', 'gene4'])
# result = pathway_pca_test(pathway_list, sc_counts)
# print(result[0])  # df
# print(result[1])  # vscore
# print(result[2])  # pa_remove


def pathway_pcascore_run(Pagwas=None, Pathway_list=None, min_pathway_size=10, max_pathway_size=1000):
    """
    计算细胞类型和单个细胞的路径PCA分数。
    """
    # 过滤路径长度
    Pa_len = {k: len(v) for k, v in Pathway_list.items()}
    Pathway_list = {k: v for k, v in Pathway_list.items() if min_pathway_size <= Pa_len[k] <= max_pathway_size}
    # 验证路径的总基因
    pa_gene = set(g for genes in Pathway_list.values() for g in genes)
    if len(set(Pagwas['avg_expr_matrix'].index).intersection(pa_gene)) < len(pa_gene) * 0.1:
        raise ValueError("Single_data和Pathway_list之间的交叉基因少于10%，请检查基因名称")
    pana = [name for name in Pathway_list if len(set(Pathway_list[name]).intersection(Pagwas['avg_expr_matrix'].index)) > 2]
    Pathway_list = {name: Pathway_list[name] for name in pana}
    Pagwas['rawPathway_list'] = Pathway_list
    # 过滤单细胞中无表达的基因
    celltypes = Pagwas['Celltype_anno']['annotation'].unique()
    data_mat = Pagwas['data_mat']
    # 假设 cellnames 是列名，annotation 是细胞类型的标签
    cellnames = Pagwas['Celltype_anno']['cellnames'].values
    annotations = Pagwas['Celltype_anno']['annotation'].values
    pana_list = []
    for celltype in celltypes:
        # 获取当前细胞类型的列索引
        celltype_indices = np.where(annotations == celltype)[0]
        # 根据索引提取数据
        scCounts = data_mat.iloc[:, celltype_indices]
        # 过滤掉和为 0 的行
        scCounts = scCounts.loc[scCounts.sum(axis=1) != 0, :]
        # 获取基因名（假设基因名是行索引）
        proper_gene_names = scCounts.index.values  # 获取基因名
        pana = [name for name in Pathway_list if len(set(Pathway_list[name]).intersection(proper_gene_names)) > 2]
        pana_list.append(pana)


    # 计算PCA分数
    print("* Start to get Pathway SVD score!")
    scPCAscore_list = []

    for celltype in tqdm(celltypes, desc='Cell Types'):
        scCounts = Pagwas['data_mat'].loc[:, Pagwas['Celltype_anno']['cellnames'][Pagwas['Celltype_anno']['annotation'] == celltype]]
        scPCAscore = pathway_pca_test(Pathway_list, scCounts)
        scPCAscore_list.append(scPCAscore)
        print(celltype)
    # 移除路径
    pa_remove = set(pa for scPCAscore in scPCAscore_list for pa in scPCAscore[2])
    # 合并所有PCA分数列表
    pca_df_list = []
    for i, scPCAscore in enumerate(scPCAscore_list):
        df = scPCAscore[0]
        if pa_remove:
            df = df[~df['name'].isin(pa_remove)]
        df['celltype'] = celltypes[i]
        pca_df_list.append(df)

    pca_df = pd.concat(pca_df_list, axis=0)
    pca_scoremat = pca_df.pivot(index='name', columns='celltype', values='score')

    pca_cell_df = pca_scoremat.fillna(0).values
    rownames_pca_cell_df = pca_scoremat.index
    colnames_pca_cell_df = pca_scoremat.columns

    pca_scCell_mat = np.concatenate([scPCAscore[1] for scPCAscore in scPCAscore_list], axis=0)
    pca_scCell_mat = pd.DataFrame(pca_scCell_mat, index=rownames_pca_cell_df, columns=Pagwas['data_mat'].columns)

    Pagwas['pca_scCell_mat'] = pca_scCell_mat
    Pagwas['merge_scexpr'] = pca_cell_df
    Pagwas['VariableFeatures'] = list(set(Pagwas['data_mat'].index).intersection(Pagwas['VariableFeatures']))

    Pagwas['pca_cell_df'] = pd.DataFrame(pca_cell_df, index=rownames_pca_cell_df, columns=colnames_pca_cell_df)

    if pa_remove:
        Pagwas['Pathway_list'] = {name: Pathway_list[name] for name in rownames_pca_cell_df}
        Pagwas['rawPathway_list'] = {name: Pathway_list[name] for name in rownames_pca_cell_df}

    return Pagwas

# 示例使用
# Pagwas = {
#     'data_mat': pd.DataFrame(np.random.rand(100, 10), index=[f'gene{i}' for i in range(100)], columns=[f'cell{i}' for i in range(10)]),
#     'VariableFeatures': [f'gene{i}' for i in range(50)],
#     'Celltype_anno': pd.DataFrame({'annotation': ['type1'] * 5 + ['type2'] * 5, 'cellnames': [f'cell{i}' for i in range(10)]})
# }
# Pathway_list = {'pathway1': [f'gene{i}' for i in range(10)], 'pathway2': [f'gene{i}' for i in range(10, 20)]}
# result = pathway_pcascore_run(Pagwas, Pathway_list)
# print(result)
