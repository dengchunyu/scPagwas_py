from setuptools import setup, find_packages

setup(
    name="scPagwas_py",
    version="0.1.0",
    packages=find_packages(where="src"),
    package_dir={"": "src"},
    install_requires=[
        # 在这里列出依赖，例如 "numpy>=1.21.0"
        "numpy","pandas","sklearn","pybedtools","dask","statsmodels","warnings","scipy","joblib"
    ],
)

