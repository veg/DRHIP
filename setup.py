from setuptools import setup, find_packages

setup(
    name="hyphy_results_toolkit",
    version="0.1.0",
    packages=find_packages(),
    install_requires=[
        "numpy>=1.20.0",
        "pandas>=1.3.0",
        "scipy>=1.7.0",
    ],
    entry_points={
        'console_scripts': [
            'hyphy-results=hyphy_results_toolkit.cli:main',
        ],
    },
    author="Danielle Callan, Hannah Verdonk, Sergei Pond",
    author_email="spond@temple.edu",  # Primary contact
    description="Tools for analyzing HyPhy results",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    url="https://github.com/veg/hyphy-results-toolkit",
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.6",
)
