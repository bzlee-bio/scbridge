"""
Setup script for scbridge Python package
"""
from setuptools import setup, find_packages
from pathlib import Path

# Read version
version = {}
with open("scbridge/_version.py") as f:
    exec(f.read(), version)

# Read README
this_directory = Path(__file__).parent
long_description = (this_directory / "README.md").read_text()

setup(
    name="scbridge",
    version=version['__version__'],
    description="Cross-platform single-cell RNA-seq data storage for Python (AnnData) and R (Seurat/SCE)",
    long_description=long_description,
    long_description_content_type="text/markdown",
    author="scbridge contributors",
    license="MIT",
    packages=find_packages(),
    install_requires=[
        "anndata>=0.8.0",
        "pandas>=1.5.0",
        "numpy>=1.21.0",
        "scipy>=1.7.0",
        "pyarrow>=10.0.0",
    ],
    python_requires=">=3.8",
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
)
