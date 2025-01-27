import pandas as pd
import numpy as np
from dask import delayed, compute
from dask.diagnostics import ProgressBar
import numpy as np
import pandas as pd
from sklearn.linear_model import LinearRegression
from statsmodels.regression.linear_model import OLS
from statsmodels.tools.tools import add_constant
import warnings
from statsmodels.stats.multitest import multipletests

def link_pathway_blocks_gwas(Pagwas, chrom_ld=None, singlecell=True, celltype=True, n_cores=1):
    # Split pathway blocks by chromosome
    Pachrom_block_list = {
        pathway: {
            chrom: group.reset_index(drop=True)  # 将每个分组转换为 DataFrame
            for chrom, group in pathway_blocks.groupby('chrom')
        }
        for pathway, pathway_blocks in Pagwas['pathway_blocks'].items()
        }


    # Split GWAS data by chromosome
    chrom_gwas_list = {
        chrom: gwas_data.sort_values(by='pos') for chrom, gwas_data in Pagwas['gwas_data'].groupby('chrom')
    }

    # Run the pathway block function with Dask for parallel processing
    with ProgressBar():
        Pagwas = pathway_block_func(
            Pagwas=Pagwas,
            Pachrom_block_list=Pachrom_block_list,
            chrom_ld=chrom_ld,
            chrom_gwas_list=chrom_gwas_list,
            singlecell=singlecell,
            celltype=celltype,
            n_cores=n_cores
        )
    
    return Pagwas


def pathway_block_func(Pagwas, Pachrom_block_list, chrom_gwas_list, singlecell=True, celltype=True, chrom_ld=None, n_cores=1):
    pathway_sclm_results = {}
    pathway_ld_gwas_data = {}

    print(f"* Start to link gwas and pathway block annotations for {len(Pachrom_block_list)} pathways!")

    # Process each pathway in parallel
    results = []
    for pathway, Pa_chrom_block in Pachrom_block_list.items():
        result = process_pathway_block(
            pathway, Pa_chrom_block, chrom_ld, chrom_gwas_list, Pagwas, singlecell, celltype, n_cores
        )
        results.append(result)
    
    # Collect results
    for pathway, pa_block, singlecell_result in results:
        if singlecell and singlecell_result is not None:
            pathway_sclm_results[pathway] = singlecell_result
        if celltype:
            pathway_ld_gwas_data[pathway] = pa_block

    # Store results in Pagwas
    if singlecell:
        Pagwas['Pathway_sclm_results'] = pd.DataFrame(pathway_sclm_results)
        Pagwas['Pathway_sclm_results'].index = Pagwas['data_mat'].columns
    
    if celltype:
        Pagwas['Pathway_ld_gwas_data'] = pathway_ld_gwas_data

    return Pagwas


def process_pathway_block(pathway, Pa_chrom_block, chrom_ld, chrom_gwas_list, Pagwas, singlecell, celltype, n_cores):
    Pa_chrom_data = []
    for chrom, chrom_block in Pa_chrom_block.items():
        ld_data = chrom_ld.get(chrom, None)
        if ld_data is None:
            print(f"Warning: {chrom} for GWAS is missing, could be a problem!")
            continue
        gwas_data = chrom_gwas_list.get(chrom, None)
        if gwas_data is None:
            print(f"Warning: {chrom} data missing, could be a problem!")
            continue
        rsids = Pagwas['snp_gene_df'].loc[Pagwas['snp_gene_df']['label'].isin(chrom_block['label']), ['rsid', 'label']]
        rsids_gwas = pd.merge(rsids, gwas_data, on='rsid', how='inner')
        if rsids_gwas.empty:
            continue
        beta_squared = rsids_gwas['beta'] ** 2
        sub_ld = ld_data[ld_data['SNP_A'].isin(rsids_gwas['rsid'])]
        Pa_chrom_data.append((rsids_gwas, beta_squared, sub_ld, chrom_block))
        
    pa_block = {
        'block_info': pd.concat([item[3] for item in Pa_chrom_data]),
        'snps': pd.concat([item[0] for item in Pa_chrom_data]),
        'y': np.concatenate([item[1] for item in Pa_chrom_data]),
        'ld_data': pd.concat([item[2] for item in Pa_chrom_data])
        }
    snp_a = set(pa_block['ld_data']['SNP_A'].dropna().values)
    snp_b = set(pa_block['ld_data']['SNP_B'].dropna().values)
    ld_data_unique = snp_a.union(snp_b)
    rsid_x = set(pa_block['snps']['rsid'].values) & ld_data_unique

    pa_block['ld_data'] = pa_block['ld_data'][pa_block['ld_data']['SNP_A'].isin(rsid_x) & pa_block['ld_data']['SNP_B'].isin(rsid_x)]
    if pa_block['ld_data'].empty:
        ld_matrix = np.eye(len(pa_block['snps']))
    else:
        ld_matrix = make_ld_matrix(pa_block['snps']['rsid'], pa_block['ld_data'])
    pa_block['ld_matrix_squared'] = ld_matrix ** 2
    pa_block['n_snps'] = len(pa_block['snps'])
    singlecell_result = None
    if singlecell:
        singlecell_result = get_pathway_sclm(pa_block, Pagwas['pca_scCell_mat'], Pagwas['data_mat'], Pagwas['rawPathway_list'],  n_cores)
    if celltype:
        pa_block = link_pwpca_block(pa_block, pca_cell_df=Pagwas['pca_cell_df'], merge_scexpr=Pagwas['avg_expr_matrix'],rawPathway_list= Pagwas['rawPathway_list'],snp_gene_df=Pagwas['snp_gene_df'])
    return pathway, pa_block, singlecell_result

def make_ld_matrix(all_snps, ld_data):
    mat_dim = len(all_snps)
    
    # 初始化对角矩阵
    ld_matrix = np.eye(mat_dim)
    
    if mat_dim == 1:
        return np.array([[1]])  # 如果只有一个 SNP，返回 1 的矩阵
    
    if mat_dim >= 2:
        # 创建一个从 all_snps 到索引的映射
        snp_index = {snp: idx for idx, snp in enumerate(all_snps)}
        
        # 设置行列名
        ld_matrix = pd.DataFrame(ld_matrix, index=all_snps, columns=all_snps)
        
        # 遍历 ld_data，填充 LD 矩阵
        for _, row in ld_data.iterrows():
            snp_a = row[0]  # SNP_A
            snp_b = row[1]  # SNP_B
            r2_value = row[2]  # r2

            # 使用 snp_index 获取对应的行列索引
            idx_a = snp_index.get(snp_a, None)
            idx_b = snp_index.get(snp_b, None)
            
            if idx_a is not None and idx_b is not None:
                ld_matrix.iloc[idx_a, idx_b] = r2_value
                ld_matrix.iloc[idx_b, idx_a] = r2_value

    return ld_matrix.values  # 返回 NumPy 数组格式的 LD 矩阵


def make_ld_matrix2(all_snps, ld_data):
    """
    Constructs an LD matrix for the specified SNPs.
    """
    snp_indices = {snp: idx for idx, snp in enumerate(all_snps)}
    mat_dim = len(all_snps)
    ld_matrix = np.eye(mat_dim)
    for _, row in ld_data.iterrows():
        i = snp_indices.get(row['SNP_A'])
        j = snp_indices.get(row['SNP_B'])
        if i is not None and j is not None:
            ld_matrix[i, j] = ld_matrix[j, i] = row['r2']
    return ld_matrix


def link_pwpca_block(pa_block, pca_cell_df, merge_scexpr, snp_gene_df, rawPathway_list):
    
    pathway = pa_block['block_info']['pathway'].unique()
    x = get_rows_by_index(pca_cell_df, pathway)

    if pa_block['snps'].shape[0] == 0:
        pa_block['include_in_inference'] = False
        pa_block['x'] = None  # Clear previous 'x'
        return pa_block[['x', 'y', 'snps', 'include_in_inference']]
    
    proper_genes = merge_scexpr.index
    genes_set = set()
    # 遍历 pathway 数组中的每个 pathway，合并基因数据
    for p in pathway:
        if p in rawPathway_list:
            # 将每个 pathway 对应的基因数据转换为 set 并合并
            genes_set.update(set(rawPathway_list[p]))

    # 如果 proper_genes 也是一个 set，取交集
    mg = list(genes_set & set(proper_genes))
    x2 = merge_scexpr.loc[mg, :]

    if x2.shape[1] == 1:
        x2 = x2.values  # Convert to matrix if there's only one column
    
    if len(mg) > 1:
        # Normalizing x2
        x2 = (x2 - x2.min(axis=0)) / (x2.max(axis=0) - x2.min(axis=0))

    if pa_block['n_snps'] > 1:
        # If more than 1 SNP
        if x2.shape[1] == 1:
            x2 = x2[pa_block['snps']['label'], :]  # Select SNPs based on label
            pa_block['n_snps'] = pa_block['snps'].shape[0]
            x = np.tile(x, (pa_block['n_snps'], 1))  # Repeat PCA scores for each SNP
            #x.index = pa_block['snps']['rsid']
            # x = x * snp_gene_df.loc[pa_block['snps']['rsid'], "slope"]
            x3 = x2 * x
        else:
            x2 = x2.loc[pa_block['snps']['label'], :]
            pa_block['n_snps'] = pa_block['snps'].shape[0]
            x = np.tile(x, (pa_block['n_snps'], 1))  # Repeat PCA scores for each SNP
            #x.index = pa_block['snps']['rsid']
            # x = x * snp_gene_df.loc[pa_block['snps']['rsid'], "slope"]
            x3 = x2 * x
    else:
        x2 = np.expand_dims(x2.loc[pa_block['snps']['label'], :], axis=0)
        x2.index = pa_block['snps']['label']
        pa_block['n_snps'] = pa_block['snps'].shape[0]
        #x.index = pa_block['snps']['rsid']
        x3 = np.multiply(x2, x)  # Element-wise multiplication

    # Removing intermediate variables
    del x
    del x2

    # Compute LD matrix squared and the final block data 'x'
    pa_block['x'] = np.dot(pa_block['ld_matrix_squared'].T, x3)
    #pa_block['x'].index = pa_block['snps']['rsid']
    #pa_block['x'].columns = merge_scexpr.columns
    # Removing intermediate variable
    del x3
    pa_block['include_in_inference'] = True
    # Returning only relevant columns
    return {key: pa_block[key] for key in ['x', 'y', 'snps', 'include_in_inference']}


def get_pathway_sclm(pa_block, pca_scCell_mat, data_mat, rawPathway_list, n_cores=1, Rns='random_name'):
    pathway = pa_block['block_info']['pathway'].unique()[0]
    
    # 获取路径索引
    row_index = np.where(pca_scCell_mat['index'] == pathway)[0]
    if row_index.size == 0:
        print(f"the '{pathway}' is not in the pca_scCell_mat")
        return pa_block 

    row_index = row_index[0]
    x = np.array(pca_scCell_mat['data'][row_index, :]).reshape(1, -1)

    if pa_block['n_snps'] == 0:
        pa_block.update({'include_in_inference': False, 'x': None})
        return pa_block

    mg = set(rawPathway_list[pathway]).intersection(set(data_mat.index))
    mg = list(mg)
    if len(mg) == 1:
        x2 = np.array(data_mat.loc[mg, :]).reshape(1, -1)
        x2 /= (x2 + 0.0001)  # Normalize
    else:
        x2 = data_mat.loc[mg, :].apply(lambda ge: ge / np.sum(ge) if np.sum(ge) > 0 else np.zeros(len(ge)), axis=0).values

    #x2 = DataMatrix(x2)
    pa_block_snps_labels = pa_block['snps']['label'].values
    # 找到与 mg 中基因匹配的行索引
    rows_to_extract = [i for i in range(len(pa_block_snps_labels)) if pa_block_snps_labels[i] in mg]

    if pa_block['n_snps'] > 1:
        x2 = x2[rows_to_extract]
        pa_block['n_snps'] = len(pa_block['snps'])
        x = np.tile(x, (pa_block['n_snps'], 1))
        x2 = x2 * x
    else:
        x2 = np.multiply(x2[rows_to_extract], x)

    pa_block['x'] = pa_block['ld_matrix_squared'] @ x2
    pa_block['include_in_inference'] = True

    noise_per_snp = pa_block['snps']['se'] ** 2

    results = None
    if pa_block['x'] is not None:
        valid_indices = (~np.isnan(pa_block['y'])) & (~np.isnan(noise_per_snp)) & \
                        ~np.isnan(pa_block['x']).any(axis=1)
        if valid_indices.any():
            results = sc_parameter_regression(
                Pagwas_x=pa_block['x'][valid_indices, :],
                Pagwas_y=pa_block['y'][valid_indices],
                n_cores=n_cores
            )
            results[np.isnan(results)] = 0
    return results


def sc_parameter_regression(Pagwas_x, Pagwas_y, n_cores=1):
    """
    Perform regression using linear regression on large-scale data with parallelization.
    """
    reg = LinearRegression(n_jobs=n_cores)
    reg.fit(Pagwas_x, Pagwas_y)
    return reg.coef_

def get_rows_by_index(pca_dict, row_names):
    """
    Given a list of row names (index labels), return the corresponding rows from the DataMatrix.

    :param pca_dict: Dictionary containing 'data', 'index', and 'columns'
    :param row_names: A list or array of row names (index labels) to fetch
    :return: A numpy array of the corresponding rows
    """
    # Ensure row_names is a list or array
    if isinstance(row_names, str):  # If a single string is passed
        row_names = [row_names]
    
    # Check that all row_names are in the index
    missing_rows = [row for row in row_names if row not in pca_dict['index']]
    if missing_rows:
        raise KeyError(f"Rows {', '.join(missing_rows)} not found in the DataMatrix.")

    # Get the indices of the requested rows
    row_indices = np.isin(pca_dict['index'], row_names)
    
    # Return the corresponding rows as a numpy array
    return pca_dict['data'][row_indices, :]


def Pagwas_perform_regression(Pathway_ld_gwas_data):
    print("** Start inference")
    
    # Vectorize Pagwas data
    vectorized_Pagwas_data = xy2vector(Pathway_ld_gwas_data)
    
    # Run regression
    lm_results = Parameter_regression(vectorized_Pagwas_data)
    
    return lm_results

def Parameter_regression(vectorized_Pagwas_data):
    lm_results = {}

    # Perform linear regression
    y = vectorized_Pagwas_data['y']
    x = vectorized_Pagwas_data['x']
    noise_per_snp = vectorized_Pagwas_data['noise_per_snp']
    
    # Adding constant for intercept
    x_with_intercept = add_constant(x)

    # OLS regression model
    model = OLS(y, x_with_intercept, weights=1/noise_per_snp).fit()
    lm_results['parameters'] = model.params
    annotation_names = ["Intercept"] + list(x.columns)
    lm_results['parameters'].index = annotation_names
    lm_results['model'] = model
    
    return lm_results

def Boot_evaluate(Pagwas, bootstrap_iters=200, part=0.5):
    print(f"* Starting bootstrap iteration for {bootstrap_iters} times")
    
    Boot_resultlist = []
    for i in range(bootstrap_iters):
        # Sample a subset of Pathway_ld_gwas_data for bootstrapping
        sampled_data = np.random.choice(Pagwas['Pathway_ld_gwas_data'], 
                                        int(len(Pagwas['Pathway_ld_gwas_data']) * part), 
                                        replace=False)
        
        # Vectorize the sampled data
        boot_results = Parameter_regression(xy2vector(sampled_data))
        boot_parameters = boot_results['parameters']
        
        # Get heritability contributions
        boot_heritability = Get_Pathway_heritability_contributions(Pagwas['pca_cell_df'], boot_parameters)
        
        Boot_resultlist.append({
            'boot_parameters': boot_parameters,
            'block_heritability': boot_heritability
        })
    
    # Compile results into DataFrame
    boot_parameters_list = [result['boot_parameters'] for result in Boot_resultlist]
    boot_parameters_df = pd.DataFrame(boot_parameters_list)
    
    Pagwas['bootstrap_results'] = Get_bootresults_df(boot_parameters_df.values, 
                                                     Pagwas['lm_results']['parameters'].index, 
                                                     Pagwas['lm_results']['parameters'])
    
    return Pagwas

def xy2vector(Pathway_ld_gwas_data):
    # Filter blocks for inclusion
    Pathway_ld_gwas_data = [block for block in Pathway_ld_gwas_data if block['include_in_inference']]
    
    y = np.concatenate([block['y'] for block in Pathway_ld_gwas_data])
    x = np.vstack([block['x'] for block in Pathway_ld_gwas_data])
    noise_per_snp = np.concatenate([block['snps']['se']**2 for block in Pathway_ld_gwas_data])
    
    # Remove NaN values
    na_elements = np.isnan(y) | np.isnan(noise_per_snp) | np.any(np.isnan(x), axis=1)
    
    return {
        'y': y[~na_elements], 
        'x': x[~na_elements], 
        'noise_per_snp': noise_per_snp[~na_elements]
    }

def Get_Pathway_heritability_contributions(pca_cell_df, parameters):
    if np.any(np.isnan(parameters)):
        warnings.warn("NA parameters found!")
        parameters[np.isnan(parameters)] = 0
    
    # Perform matrix multiplication to get pathway heritability contributions
    pathway_block_info = np.dot(pca_cell_df.values, parameters[pca_cell_df.columns])
    return pd.Series(pathway_block_info, index=pca_cell_df.index)

def Get_bootresults_df(value_collection, annotations, model_estimates):
    # Prepare DataFrame to summarize bootstrap results
    value_collection = np.array(value_collection) if len(value_collection.shape) == 1 else value_collection
    bootstrap_estimate = np.mean(value_collection, axis=0)
    bootstrap_error = np.std(value_collection, axis=0)
    
    bt_value = bootstrap_estimate / bootstrap_error
    bp_value = 1 - stats.norm.cdf(bt_value)  # p-value calculation
    
    # Bias-corrected estimate
    bias_corrected_estimate = 2 * model_estimates - bootstrap_estimate
    
    # Confidence intervals
    CI_lo = np.percentile(value_collection, 2.5, axis=0)
    CI_hi = np.percentile(value_collection, 97.5, axis=0)
    
    result_df = pd.DataFrame({
        'annotation': annotations,
        'bootstrap_estimate': bootstrap_estimate,
        'bootstrap_error': bootstrap_error,
        'bt_value': bt_value,
        'bp_value': bp_value,
        'bias_corrected_estimate': bias_corrected_estimate,
        'CI_lo': CI_lo,
        'CI_hi': CI_hi
    })
    
    return result_df


def celltype_regression(Pagwas, iters_celltype, output_dirs, output_prefix):
    """
    Perform regression and bootstrapping for pathway block data, then adjust p-values using FDR.

    Parameters:
    Pagwas (dict): A dictionary containing the necessary input data such as `Pathway_ld_gwas_data`, etc.
    celltype (bool): A flag to indicate whether to perform the celltype-specific regression.
    iters_celltype (int): Number of bootstrap iterations.
    output_dirs (str): Directory to save the output CSV file.
    output_prefix (str): Prefix for the output CSV file name.

    Returns:
    dict: The updated Pagwas dictionary with regression and bootstrap results.
    """
    # Perform regression
    Pagwas['lm_results'] = Pagwas_perform_regression(
        Pathway_ld_gwas_data=Pagwas['Pathway_ld_gwas_data']
    )
    
    # Perform bootstrap evaluation
    Pagwas = Boot_evaluate(Pagwas, bootstrap_iters=iters_celltype, part=0.5)

    # Adjust p-values using FDR (Benjamini-Hochberg correction)
    Pagwas['bootstrap_results']['bp_value'] = multipletests(
        Pagwas['bootstrap_results']['bp_value'], method='fdr_bh')[1]

    # Clean up unused data
    Pagwas['Pathway_ld_gwas_data'] = None
    Pagwas['merge_scexpr'] = None

    # Save the bootstrap results to a CSV file
    bootstrap_results = Pagwas['bootstrap_results']
    output_file = f"./{output_dirs}/{output_prefix}_celltypes_bootstrap_results.csv"
    bootstrap_results.to_csv(output_file, index=False, quoting=2)  # quote=False for CSV
    
    return Pagwas
