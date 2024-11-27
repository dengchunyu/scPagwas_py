from setuptools import setup, find_packages

setup(
    name="scPagwas_py",
    version="0.1.0",
    packages=find_packages(where="src"),
    package_dir={"": "src"},
    install_requires=[
    "numpy",
    "pandas",
    "scikit-learn",
    "pybedtools",
    "dask",
    "statsmodels",
    "scipy",
    "joblib",
    # warnings 不是一个独立的包，删除它
    ],
)

