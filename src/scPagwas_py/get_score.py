import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix
from sklearn.linear_model import LinearRegression
from scipy.stats import wilcoxon
from statsmodels.stats.multitest import multipletests
from joblib import Parallel, delayed
from scipy.stats import mannwhitneyu
from scipy.stats import norm
from scipy import stats
from statsmodels.stats import multitest

def get_pathway_sclm(pa_block, pca_scCell_mat, data_mat, rawPathway_list, n_cores=1, backingpath='', Rns=''):
    pathway = np.unique(pa_block['block_info']['pathway'])
    row_indices = [np.where(pca_scCell_mat['index'] == name)[0] for name in pathway]
    row_indices = [index[0] for index in row_indices if index.size > 0]  # 只保留存在的行索引

    x = np.array(pca_scCell_mat['data'][row_indices, :]).reshape(1, -1)

    if pa_block['n_snps'] == 0:
        pa_block['include_in_inference'] = False
        pa_block['x'] = None
        return pa_block

    mg = np.intersect1d(rawPathway_list[pathway], data_mat.index)
    if len(mg) == 1:
        x2 = data_mat.loc[mg].values.reshape(1, -1)
        x2 /= (x2 + 0.0001)
    else:
        x2 = data_mat.loc[mg].apply(lambda ge: ge / ge.sum() if ge.sum() > 0 else np.zeros(len(ge)), axis=0).values

    x2 = csr_matrix(x2)

    if pa_block['n_snps'] > 1:
        x2 = x2[pa_block['snps']['label'], :]
        pa_block['n_snps'] = x2.shape[0]
        x = np.repeat(x, pa_block['n_snps'], axis=0)

        x2 = x2.multiply(x)
    else:
        x2 = x2[pa_block['snps']['label'], :].reshape(1, -1)
        pa_block['n_snps'] = x2.shape[0]
        x = np.repeat(x, pa_block['n_snps'], axis=0)
        x2 = x2.multiply(x)

    pa_block['x'] = csr_matrix(pa_block['ld_matrix_squared'].dot(x2))
    pa_block['include_in_inference'] = True

    noise_per_snp = pa_block['snps']['se'] ** 2

    if pa_block['x'] is not None:
        if pa_block['n_snps'] > 2:
            na_elements = np.isnan(pa_block['y']) | np.any(np.isnan(pa_block['x'].toarray()), axis=1) | np.isnan(noise_per_snp)

            results = sc_parameter_regression(
                Pagwas_x=pa_block['x'][~na_elements],
                Pagwas_y=pa_block['y'][~na_elements],
                noise_per_snp=noise_per_snp[~na_elements],
                Rns=Rns,
                backingpath=backingpath,
                n_cores=n_cores
            )
            results[np.isnan(results)] = 0
            return pd.Series(results, index=data_mat.columns)
        else:
            return None
    else:
        return None


def sc_pagwas_perform_score(Pagwas, Single_data, iters_singlecell,n_topgenes,remove_outlier=True,random_times=100,n_jobs=1):
    Pathway_sclm_results = Pagwas['Pathway_sclm_results']
    Pathway_names = Pathway_sclm_results.columns

    pathway_expr = {}
    for pa in Pathway_names:
        a = np.intersect1d(Pagwas['Pathway_list'][pa], Pagwas['data_mat'].index)
        if len(a) == 0:
            continue
        elif len(a) == 1:
            pathway_expr[pa] = Pagwas['data_mat'].loc[a].values
        else:
            pathway_expr[pa] = Pagwas['data_mat'].loc[a].mean(axis=0).values

    pathway_expr_df = pd.DataFrame.from_dict(pathway_expr, orient='index').T
    row_indices = [np.where(Pagwas['pca_scCell_mat']['index'] == name)[0] for name in Pathway_names]
    row_indices = [index[0] for index in row_indices if index.size > 0]  # 只保留存在的行索引

    pa_exp_mat = Pagwas['pca_scCell_mat']['data'][row_indices, :].T * pathway_expr_df.values
    Pagwas['Pathway_single_results'] = Pathway_sclm_results[Pathway_names].values * pa_exp_mat

    #cl = np.unique(Pagwas['Celltype_anno']['annotation'])
    scPagwas_gPAS_score = Pagwas['Pathway_single_results'].sum(axis=1)

    if remove_outlier:
        scPagwas_gPAS_score = sc_pagwas_score_filter(scPagwas_gPAS_score)

    Pagwas['scPagwas.gPAS.score'] = scPagwas_gPAS_score
    Pagwas = sc_get_pcc2(Pagwas, random_times=random_times,n_jobs=n_jobs)
    Pagwas = get_trs_score(Pagwas, Single_data, iters_singlecell=iters_singlecell, n_topgenes=n_topgenes)

    return Pagwas


def sc_parameter_regression(Pagwas_x, Pagwas_y):
    if Pagwas_x.shape[1] <= 20000:
        model = LinearRegression()
        model.fit(Pagwas_x, Pagwas_y)
        return model.coef_
    else:
        # Split Pagwas_x into chunks of 10,000 columns
        n = Pagwas_x.shape[1] // 10000
        split_cols = np.array_split(np.arange(Pagwas_x.shape[1]), n)
        results = [LinearRegression().fit(Pagwas_x[:, cols], Pagwas_y).coef_ for cols in split_cols]
        return np.concatenate(results)


def sc_pagwas_score_filter(scPagwas_score):
    scPagwas_score = np.nan_to_num(scPagwas_score, nan=0)
    lower_bound, upper_bound = np.percentile(scPagwas_score, [1, 99])
    scPagwas_score[scPagwas_score < lower_bound] = lower_bound
    scPagwas_score[scPagwas_score > upper_bound] = upper_bound
    return scPagwas_score


def pa_meanexpr(Pagwas):
    Pathway_names = Pagwas['Pathway_list'].keys()
    pathway_expr = {}
    for pa in Pathway_names:
        a = np.intersect1d(Pagwas['Pathway_list'][pa], Pagwas['data_mat'].index)
        if len(a) == 0:
            continue
        elif len(a) == 1:
            pathway_expr[pa] = Pagwas['data_mat'].loc[a].values
        else:
            pathway_expr[pa] = Pagwas['data_mat'].loc[a].mean(axis=0).values

    pathway_expr_df = pd.DataFrame.from_dict(pathway_expr, orient='index').T

    row_indices = [np.where(Pagwas['pca_scCell_mat']['index'] == name)[0] for name in Pathway_names]
    row_indices = [index[0] for index in row_indices if index.size > 0]  # 只保留存在的行索引
    pa_exp_mat = Pagwas['pca_scCell_mat']['data'][row_indices, :].T.values * pathway_expr_df.values

    Pagwas['pa_exp_mat'] = pa_exp_mat
    return Pagwas


def merge_gPas(Pagwas, pmat_merge):
    if 'pa_exp_mat' not in Pagwas:
        Pagwas = pa_meanexpr(Pagwas)
        
    pa_exp_mat = Pagwas['pca_scCell_mat'].data.T * Pagwas['pa_exp_mat']
    Pathway_names = np.intersect1d(pmat_merge.columns, Pagwas['pa_exp_mat'].columns)
    Pathway_single_results = pmat_merge[Pathway_names] * pa_exp_mat[Pathway_names]
    
    scPagwas_gPAS_score = Pathway_single_results.sum(axis=1)
    scPagwas_gPAS_score = sc_pagwas_score_filter(scPagwas_gPAS_score)
    return scPagwas_gPAS_score




def cor_sparse(X, Y):
    """
    Compute the correlation matrix between two matrices X and Y using a sparse correlation method.
    
    Parameters:
    X (numpy.ndarray): A data matrix with shape (n_samples, n_features).
    Y (numpy.ndarray): Another data matrix with shape (n_samples, n_features).
    
    Returns:
    numpy.ndarray: The correlation matrix between X and Y.
    """
    # Ensure X and Y have the same number of rows
    assert X.shape[0] == Y.shape[0], "X and Y must have the same number of rows (samples)."
    
    # Number of samples (rows)
    n = X.shape[0]
    
    # Mean of columns in X and Y
    muX = np.mean(X, axis=0)
    muY = np.mean(Y, axis=0)
    
    # Covariance matrix calculation
    covmat = (X.T @ Y - n * np.outer(muX, muY))
    
    # Standard deviations of X and Y
    sdvecX = np.sqrt(np.sum(X**2, axis=0) - n * muX**2)
    sdvecY = np.sqrt(np.sum(Y**2, axis=0) - n * muY**2)
    
    # Compute correlation matrix (element-wise division)
    cormat = covmat / np.outer(sdvecX, sdvecY)
    
    return cormat


def compute_random_correlation(datamat, gPas, select_num, random_times=100,n_jobs=1):
    sparse_cor_list = []

    # 预先生成所有的随机索引
    random_indices = [np.random.choice(datamat.shape[1], select_num, replace=False) for _ in range(random_times)]
    def cor_sparse(X, Y):
        assert X.shape[0] == Y.shape[0], "X and Y must have the same number of rows (samples)."
        n = X.shape[0]  
        muX = np.mean(X, axis=0)
        muY = np.mean(Y, axis=0)
        covmat = (X.T @ Y - n * np.outer(muX, muY))
        sdvecX = np.sqrt(np.sum(X**2, axis=0) - n * muX**2)
        sdvecY = np.sqrt(np.sum(Y**2, axis=0) - n * muY**2)
        cormat = covmat / np.outer(sdvecX, sdvecY)
        return cormat
    # 并行计算相关性
    def compute_single_correlation(index):
        gPas_select = gPas[index]
        datamat_selected = datamat.iloc[:, index].to_numpy().T
        sparse_cor = cor_sparse(datamat_selected, gPas_select)
        sparse_cor = np.nan_to_num(sparse_cor, nan=0.0)  # Clean up NaN and Inf
        return sparse_cor[:, 0]
    
    # 并行处理每次采样
    sparse_cor_list = Parallel(n_jobs=n_jobs)(delayed(compute_single_correlation)(index) for index in random_indices)
    sparse_cor_df = pd.DataFrame(sparse_cor_list).T
    sparse_cor_mean = np.mean(sparse_cor_df, axis=1)
    return sparse_cor_mean

def sc_get_pcc2(Pagwas, random_times=100,n_jobs=1):
    """
    Compute PCC for each gene pathway in Pagwas object and adjust p-values using Wilcoxon test.

    Parameters:
    Pagwas (dict): Pagwas result dictionary, contains 'scPagwas.gPAS.score' and 'data_mat'.
    random_times (int): Number of random iterations for PCC computation.

    Returns:
    dict: Updated Pagwas dictionary with PCC and p-value information.
    """
    scPagwas_gPAS_score = Pagwas['scPagwas.gPAS.score']
    data_mat = Pagwas['data_mat']
    pcc = compute_random_correlation(gPas=scPagwas_gPAS_score, datamat=data_mat,
                     random_times=random_times, select_num=int(data_mat.shape[1] / 5),n_jobs=n_jobs)
    
    # Create DataFrame for PCC
    pcc_df = pd.DataFrame(pcc, columns=['PCC'])
    
    # Identify top 5% cells based on gPAS score
    top5_idx = np.argsort(scPagwas_gPAS_score)[::-1][:int(0.05 * len(scPagwas_gPAS_score))]
    #top5_cellnames = np.array(list(scPagwas_gPAS_score.keys()))[top5_idx]

    # Perform Wilcoxon test for each gene and adjust p-value using Bonferroni correction
    # 确保 data_mat 是一个 numpy 数组
    gene_names=data_mat.index
    if isinstance(data_mat, pd.DataFrame):
        data_mat = data_mat.to_numpy()

    # 定义 p_list 用于存储每个基因的 p 值
    p_list = []

    # 进行 Wilcoxon 测试
    for gene_idx in range(data_mat.shape[0]):
        group_top5 = data_mat[gene_idx, top5_idx]
        group_rest = data_mat[gene_idx, ~np.isin(np.arange(data_mat.shape[1]), top5_idx)]
        # 执行 Wilcoxon 检验
        _, p_value = mannwhitneyu(group_top5, group_rest, alternative='greater')
        p_list.append(p_value)
    # Adjust p-values using Bonferroni correction
    p_list_adjusted = multipletests(p_list, method='bonferroni')[1]

    # Add p-values and adjusted p-values to DataFrame
    pcc_df['pvalue'] = p_list
    pcc_df['adj_pvalue'] = p_list_adjusted
    pcc_df['adj_logp'] = -np.log10(pcc_df['adj_pvalue'])

    # Handle any invalid logp values (e.g., Inf or NaN)
    pcc_df['adj_logp'] = np.where(np.isfinite(pcc_df['adj_logp']), pcc_df['adj_logp'], 
                                  pcc_df['adj_logp'].max() + 1)

    # Compute weighted PCC
    pcc_df['weight_pcc'] = pcc_df['adj_logp'] * pcc_df['PCC']

    # Store the results in Pagwas dictionary
    Pagwas['PCC'] = pcc_df
    Pagwas['PCC'].index=gene_names

    return Pagwas


def get_correct_bg_p(data_mat, scPagwas_TRS_Score, iters_singlecell, n_topgenes, scPagwas_topgenes):
    # 获取 `data_mat` 和 `scPagwas_topgenes` 的交集基因
    genes = np.intersect1d(data_mat.index, scPagwas_topgenes)
    if len(genes) == 0:
        raise ValueError("No common genes found between `data_mat` and `scPagwas_topgenes`.")
    
    # 初始化控制评分矩阵
    gene_matrix = data_mat.to_numpy()  # 转换为 NumPy 数组
    mat_ctrl_raw_score = np.zeros((gene_matrix.shape[1], iters_singlecell))
    dic_ctrl_list = []
    
    # 在迭代中随机选择控制基因集并计算控制评分
    for i in range(iters_singlecell):
        np.random.seed(i)
        # 从 `data_mat.index` 中随机选择控制基因
        selected_genes = np.random.choice(data_mat.index, n_topgenes, replace=False)
        dic_ctrl_list.append(selected_genes)
        # 获取控制基因在 `data_mat` 中的索引
        control_gene_indices = data_mat.index.get_indexer(selected_genes)
        # 在 `gene_matrix` 中计算选中控制基因的均值评分
        mat_ctrl_raw_score[:, i] = gene_matrix[control_gene_indices, :].mean(axis=0)
    
    # 计算基因方差
    df_gene_var = pd.DataFrame({'gene': data_mat.index, 'var': gene_matrix.var(axis=1)}).set_index('gene')
    var_ratio_c2t = np.array([
        df_gene_var.loc[ctrl_genes, 'var'].sum() / df_gene_var.loc[genes, 'var'].sum()
        for ctrl_genes in dic_ctrl_list
    ])
    
    # 计算 p 值和调整后的分布
    correct_pdf = correct_background(scPagwas_TRS_Score, mat_ctrl_raw_score, var_ratio_c2t)
    correct_pdf.index = data_mat.columns
    return correct_pdf


def correct_background(scPagwas_TRS_Score, mat_ctrl_raw_score, var_ratio_c2t):
    scPagwas_TRS_Score = scPagwas_TRS_Score - np.mean(scPagwas_TRS_Score)
    mat_ctrl_raw_score = (mat_ctrl_raw_score - np.mean(mat_ctrl_raw_score, axis=0)) / np.sqrt(var_ratio_c2t)
    
    mean_ctrl = np.mean(mat_ctrl_raw_score, axis=1)
    std_ctrl = np.std(mat_ctrl_raw_score, axis=1)
    
    norm_scores = (scPagwas_TRS_Score - mean_ctrl) / std_ctrl
    mat_ctrl_norm_scores = (mat_ctrl_raw_score - mean_ctrl[:, None]) / std_ctrl[:, None]
    
    min_norm_score = min(np.min(norm_scores), np.min(mat_ctrl_norm_scores))
    norm_scores[scPagwas_TRS_Score == 0] = min_norm_score - 1e-3
    mat_ctrl_norm_scores[mat_ctrl_raw_score == 0] = min_norm_score
    
    pooled_p_values = get_p_from_empi_null(norm_scores, mat_ctrl_norm_scores.flatten())
    
    # 使用 statsmodels 进行多重检验校正
    _, adj_p_values, _, _ = multitest.multipletests(pooled_p_values, method='fdr_bh')
    
    correct_pdf = pd.DataFrame({
        'raw_score': scPagwas_TRS_Score,
        'pooled_p': pooled_p_values,
        'adj_p': adj_p_values,
        'nlog10_pooled_p': -np.log10(pooled_p_values),
        'pooled_z': np.clip(-stats.norm.ppf(pooled_p_values), -10, 10)
    })
    
    return correct_pdf

def get_p_from_empi_null(v_t, v_t_null):
    v_pos = np.sum(v_t_null <= v_t[:, None], axis=1)
    v_p = (len(v_t_null) - v_pos + 1) / (len(v_t_null) + 1)
    return v_p

def merge_celltype_p(single_p, celltype):
    celltype_p = pd.DataFrame({'celltype': celltype, 'pvalue': single_p})
    # 显式设置 observed=False 以消除警告
    celltype_p = celltype_p.groupby('celltype', observed=False)['pvalue'].apply(p_merge).reset_index()
    return celltype_p

def p_merge(pvalues):
    zvalues = -np.sqrt(2) * stats.norm.ppf(pvalues / 2)
    ztotal = np.mean(zvalues)
    p_total = stats.norm.cdf(-np.abs(ztotal))
    return p_total


def get_trs_score(Pagwas, Single_data, iters_singlecell, n_topgenes):
    """
    Perform the full scPagwas analysis pipeline.
    
    Parameters:
    - Pagwas: Dictionary containing the PCC and other necessary data
    - Single_data: DataFrame, Single cell data (Assumed to have the same structure as Seurat object)
    - iters_singlecell: Number of iterations for the background correction
    - n_topgenes: Number of top genes to consider for the analysis
    
    Returns:
    - Pagwas: Updated dictionary with the results
    """

    # 根据 'PCC' 列排序，并选择排序后前 n_topgenes 的基因索引
    scPagwas_topgenes = Pagwas['PCC'].sort_values(by='PCC', ascending=False).head(n_topgenes).index

    # 2. AddModuleScore 
    # Create a module score by averaging the expression levels of top genes and down-regulated genes
    TRS_Score = add_module_score(data_mat=Pagwas['data_mat'], features=scPagwas_topgenes, nbin=24, ctrl=100)

    # 3. Correct background using Get_CorrectBg_p function
    correct_pdf = get_correct_bg_p(
        data_mat=Pagwas['data_mat'],
        scPagwas_TRS_Score=TRS_Score,
        iters_singlecell=iters_singlecell,
        n_topgenes=n_topgenes,
        scPagwas_topgenes=scPagwas_topgenes
    )
    Pagwas['Random_Correct_BG_pdf'] = correct_pdf

    # 4. Merge p-values by celltype
    print("* Get Merged pvalue for each celltype!")
    Pagwas['Merged_celltype_pvalue'] = merge_celltype_p(single_p=correct_pdf['pooled_p'], celltype=Pagwas['Celltype_anno']['annotation'])

    # 5. Update the scPagwas.TRS.Score
    Pagwas['scPagwas.TRS.Score'] = TRS_Score

    # 6. Prepare final data for output
    result_df = pd.DataFrame({
        'scPagwas.TRS.Score': Pagwas['scPagwas.TRS.Score'],
        'scPagwas.gPAS.score': Pagwas['scPagwas.gPAS.score'],
        'Random_Correct_BG_p': correct_pdf['pooled_p'],
        'Random_Correct_BG_adjp': correct_pdf['adj_p'],
        'Random_Correct_BG_z': correct_pdf['pooled_z']
    })
    # Save the merged celltype p-values
    merged_celltype_pvalue = Pagwas['Merged_celltype_pvalue']
    merged_celltype_pvalue.to_csv("Merged_celltype_pvalue.csv", index=False)
    # Optionally save the result dataframe
    result_df.to_csv("cell_trs_gpas_p_result.csv", index=False)

    return Pagwas


def add_module_score(data_mat, features, pool=None, nbin=24, ctrl=100, seed=1):
    np.random.seed(seed)
    
    # 使用数据框的数组表示，行是基因，列是细胞
    assay_data = data_mat.to_numpy()
    pool = pool or data_mat.index.tolist()
    
    # 筛选基因池并计算每个基因的均值
    pool_data = assay_data[data_mat.index.isin(pool), :]
    data_avg = np.mean(pool_data, axis=1)  # 按行计算均值，即每个基因的均值
    data_cut = pd.cut(data_avg + np.random.normal(scale=1e-30, size=data_avg.shape), bins=nbin, labels=False, right=False)
    data_cut = pd.Series(data_cut, index=data_mat.index[data_mat.index.isin(pool)])  # 基因池的基因名作为索引
    
    # 获取控制基因
    ctrl_genes = []
    for gene in features:
        if gene in data_cut.index:
            # 获取与当前基因相同分组的基因
            similar_genes = data_cut[data_cut == data_cut[gene]].index
            if len(similar_genes) > 0:
                # 随机选择控制基因
                ctrl_genes.extend(np.random.choice(similar_genes, size=min(ctrl, len(similar_genes)), replace=False))
    
    ctrl_genes = list(set(ctrl_genes))  # 去重处理
    if not ctrl_genes:
        raise ValueError("No control genes selected. Check if `features` match `data_mat.index` or if data_cut has valid similar genes.")
    
    # 计算控制基因评分和特征评分
    ctrl_scores = np.mean(assay_data[data_mat.index.isin(ctrl_genes), :], axis=0)
    feature_scores = np.mean(assay_data[data_mat.index.isin(features), :], axis=0)
    
    # 计算差值评分
    feature_scores_use = feature_scores - ctrl_scores    
    return feature_scores_use
