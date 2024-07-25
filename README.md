# scPagwas_py
## Data Input

This is a simple example package. You can use
[GitHub-flavored Markdown](https://guides.github.com/features/mastering-markdown/)
to write your content.

### single cell data
```python
import scanpy as sc
import pandas as pd

mtx_file_counts="/Users/chunyudeng/Library/CloudStorage/OneDrive-共享的库-Onedrive/RPakage/scPagwas_py/data/GSE115978_counts.csv.gz"
anno_file="/Users/chunyudeng/Library/CloudStorage/OneDrive-共享的库-Onedrive/RPakage/scPagwas_py/data/GSE115978_cell.annotations.csv.gz"
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
adata.write("/Users/chunyudeng/Library/CloudStorage/OneDrive-共享的库-Onedrive/RPakage/scPagwas_py/data/scdata.adata")
```

### other data

```python
#gwas data
import pandas as pd
import os
import csv
os.chdir('/Users/chunyudeng/Library/CloudStorage/OneDrive-共享的库-Onedrive/RPakage/scPagwas_py/data/')


# 打开文件并指定制表符为分隔符
gwas_data = pd.read_csv('GWAS_summ_example.txt', sep=' ')

#pathway data
pathway_data = {}
with open('Genes_by_pathway_kegg.csv', mode='r', encoding='utf-8') as file:
    csv_reader = csv.reader(file)
    for row in csv_reader:
        # 第一列是路径，第二列是基因列表
        pathway = row[0]
        genes = row[1].split(',')
        pathway_data[pathway] = genes

# 打印结果
print(pathway_data)

#block_annotation

block_annotation = pd.read_csv("block_annotation.csv")

#chrom_LD
# 定义存储数据的字典
chrom_LD = {}
os.chdir('/Users/chunyudeng/Library/CloudStorage/OneDrive-共享的库-Onedrive/RPakage/scPagwas_py/data/LD')

# 读取每个CSV文件并存储到字典中
for i in range(1, 23):
    file_name = f'{i}.csv'
    chr_key = f'chr{i}'
    df = pd.read_csv(file_name)
    chrom_LD[chr_key] = df

```

## Run scPagwas

```python
adata = sc.read_h5ad("/Users/chunyudeng/Library/CloudStorage/OneDrive-共享的库-Onedrive/RPakage/scPagwas_py/data/scdata.adata")
Pagwas = scdata_process(Pagwas=None,adata=adata, cell_type_col='cell.types')
Pagwas = GWAS_summary_input(Pagwas=Pagwas, gwas_data=gwas_data)

Pathway_list=pathway_data
min_pathway_size=10
max_pathway_size=1000

```