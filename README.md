# scPagwas_py

## install the scPagwas_py

```shell
cd /share/pub/dengcy/
git clone https://github.com/dengchunyu/scPagwas_py.git
conda create -n scPagwas_env python=3.8
conda activate scPagwas_env
cd scPagwas_py
pip install .

#安装依赖包
conda install -c bioconda bedtools
pip install scanpy numpy pandas scikit-learn pybedtools dask statsmodels scipy joblib

#update packages
pip install --upgrade git+https://github.com/dengchunyu/scPagwas_py.git
```

## Data Input
### single cell data

```python
import scanpy as sc
import pandas as pd

mtx_file_counts="/share/pub/dengcy/scPagwas_py/data/GSE115978_counts.csv.gz"
anno_file="/share/pub/dengcy/scPagwas_py/data/GSE115978_cell.annotations.csv.gz"
adata = sc.read_csv(mtx_file_counts, first_column_names=True)
adata = adata.T
adata.var_names_make_unique()
# 读取细胞注释数据
adata_annotations = pd.read_csv(anno_file)
adata_annotations.index = adata_annotations['cells']
adata.obs = adata_annotations
#过滤
sc.pp.filter_genes(adata, min_cells=50)
sc.pp.filter_cells(adata, min_genes=200)

#标准化
sc.pp.normalize_total(adata)
sc.pp.log1p(adata)
#由于我们下载的数据自带分类标签信息，因此我们就不进行细胞聚类和注释了。
adata.write("/share/pub/dengcy/scPagwas_py/data/scdata.adata")
```

### other data

```python
#gwas data
import pandas as pd
import os
import csv

os.chdir('/share/pub/dengcy/scPagwas_py/data/')
os.chdir('/Users/chunyudeng/Library/CloudStorage/OneDrive-共享的库-Onedrive/RPakage/scPagwas_py/data/')

# 打开文件并指定制表符为分隔符
gwas_data = pd.read_csv('GWAS_summ_example.txt', sep=' ')

gwas_data = pd.read_csv('/share/pub/zhouyj/brain_TF/GWAS/AD_pagwas.txt', sep=' ')

#pathway data

pathway_data = {}
with open('Genes_by_pathway_kegg.csv', mode='r', encoding='utf-8') as file:
    csv_reader = csv.reader(file)
    for row in csv_reader:
        pathway = row[0]
        genes = row[1].split(',')
        pathway_data[pathway] = genes
#block_annotation
block_annotation = pd.read_csv("block_annotation.csv")
#chrom_LD
# 定义存储数据的字典
chrom_LD = {}
os.chdir('/share/pub/dengcy/scPagwas_py/data/LD')
os.chdir('/Users/chunyudeng/Library/CloudStorage/OneDrive-共享的库-Onedrive/RPakage/scPagwas_py/data/LD')

# 读取每个CSV文件并存储到字典中
for i in range(1, 23):
    file_name = f'{i}.csv'
    chr_key = f'chr{i}'
    df = pd.read_csv(file_name)
    chrom_LD[chr_key] = df

```

## Run the input data

```python          
#os.chdir('/share/pub/dengcy/scPagwas_py/src/scPagwas_py')
#import input_pp
#importlib.reload(input_pp)

#import regress_method
from scPagwas_py import input_pp,regress_method,get_score
#去除scPagwas_py，重新导入
import importlib
importlib.reload(input_pp)
importlib.reload(regress_method)
importlib.reload(get_score)

os.chdir("/share/pub/dengcy/scPagwas_py/data")
os.chdir("/Users/chunyudeng/Library/CloudStorage/OneDrive-共享的库-Onedrive/RPakage/scPagwas_py/data")

adata = sc.read_h5ad("scdata.adata")
Pagwas = input_pp.scdata_process(Pagwas=None,adata=adata, cell_type_col='cell.types')
Pagwas = input_pp.GWAS_summary_input(Pagwas=Pagwas, gwas_data=gwas_data,block_annotation=block_annotation)
Pagwas = input_pp.pathway_pcascore_run(Pagwas=Pagwas, Pathway_list=pathway_data, min_pathway_size=10, max_pathway_size=1000,n_jobs=10)
Pagwas = input_pp.pathway_annotation_input(Pagwas=Pagwas, block_annotation=block_annotation)
```

## Run the regression method

```python
#import importlib
#import regress_method
Pagwas = regress_method.link_pathway_blocks_gwas(Pagwas=Pagwas, chrom_ld=chrom_LD, singlecell=True, celltype=True,n_cores=5)
print("Pagwas 包含的键:", list(Pagwas.keys()))
if 'Pathway_sclm_results' in Pagwas and Pagwas['Pathway_sclm_results'] is not None:
    print("通路单细胞回归结果前 5 行:\n", Pagwas['Pathway_sclm_results'].head())
else:
    print("未找到通路单细胞回归结果。")
if 'Pathway_ld_gwas_data' in Pagwas:
    for pathway, data in Pagwas['Pathway_ld_gwas_data'].items():
        print(f"\n通路名称: {pathway}")
        print("SNP 列表:\n", data['snps'].head())
else:
    print("未找到通路 LD-GWAS 数据块。")
if 'lm_results' in Pagwas:
    print("回归参数:\n", Pagwas['lm_results']['parameters'])
    print("模型摘要:\n", Pagwas['lm_results']['model'].summary())
else:
    print("未找到线性回归结果。")
#import get_score
#importlib.reload(get_score)
Pagwas = get_score.sc_pagwas_perform_score(Pagwas=Pagwas, remove_outlier=True,Single_data=adata, random_times=100, iters_singlecell=100,n_topgenes=1000,n_jobs=5)
```