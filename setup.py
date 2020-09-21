# @Created Date: 2019-11-24 09:07:03 pm
# @Filename: setup.py
# @Email:  1730416009@stu.suda.edu.cn
# @Author: ZeFeng Zhu
# @Last Modified: 2019-12-23 04:27:14 pm
# @Copyright (c) 2019 MinghuiGroup, Soochow University
from setuptools import setup, find_packages, find_namespace_packages
import pdb_profiling

with open("README.md", "rt") as f:
    readme = f.read()


setup(
    name="pdb_profiling",
    version=pdb_profiling.__version__,

    packages=find_namespace_packages(),
    install_requires=[
          'aiosqlite',
          'aioftp',
          'aiohttp',
          'aiofiles',
          'unsync',
          'tenacity',
          'furl',
          'orjson',
          'pyexcel',
          'tablib',
          'pandas',
          'numpy',
          'textdistance',
          'databases',
          'neo4j',
          'tqdm',
          'orm',
          'smart_open'
     ],
    license="MIT",
    author_email="1730416009@stu.suda.edu.cn",
    maintainer="ZeFeng Zhu",
    maintainer_email="1730416009@stu.suda.edu.cn",
    description="Profiling Protein Structures from Protein Data Bank and integrate various resources.",
    long_description=readme,
    long_description_content_type="text/markdown",
    url="https://github.com/NatureGeorge/pdb-profiling",
    python_requires=">=3.6.*",
    classifiers=[
        "Operating System :: OS Independent",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8"
    ],
    )

"""
Packaged By
python setup.py sdist bdist_wheel
twine upload --repository-url https://test.pypi.org/legacy/ dist/*
twine upload --repository-url https://upload.pypi.org/legacy/ dist/*

> https://packaging.python.org/tutorials/packaging-projects/
"""
