import os
import sys

from setuptools import find_packages, setup

# Add the project directory to the path so we can import the version
sys.path.insert(0, os.path.abspath(os.path.dirname(__file__)))
from drhip.version import __version__

# Read the contents of README.md
this_directory = os.path.abspath(os.path.dirname(__file__))
with open(os.path.join(this_directory, "README.md"), encoding="utf-8") as f:
    long_description = f.read()

setup(
    name="drhip",  # DRHIP: Data Reduction for HyPhy with Inference Processing
    version=__version__,
    packages=find_packages(include=["drhip", "drhip.*"]),
    install_requires=[
        "numpy>=1.20.0",
        "pandas>=1.3.0",
        "scipy>=1.7.0",
    ],
    extras_require={
        "dev": [
            "pytest>=7.0.0",
            "pytest-cov>=4.0.0",
            "pre-commit>=3.0.0",
            "black>=24.0.0",
            "isort>=5.12.0",
            "ruff>=0.1.0",
        ],
    },
    entry_points={
        "console_scripts": [
            "drhip=drhip.cli:main",
        ],
    },
    author="Danielle Callan, Hannah Verdonk, Sergei L Kosakovsky Pond",
    author_email="spond@temple.edu",
    description="DRHIP: Data Reduction for HyPhy with Inference Processing - A toolkit for analyzing and summarizing HyPhy evolutionary selection analysis results",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/veg/drhip",
    project_urls={
        "Bug Reports": "https://github.com/veg/drhip/issues",
        "Source": "https://github.com/veg/drhip",
    },
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.7",
    include_package_data=True,
    zip_safe=False,
)
