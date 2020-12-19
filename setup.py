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
    entry_points={'console_scripts': ['pdb_profiling=pdb_profiling.commands.command:Interface']},
    install_requires=[
        'aiosqlite>=0.13.0',
        'aiohttp>=3.6.2',
        'aiofiles>=0.6.0',
        'unsync>=1.2.1',
        'tenacity>=6.3.0',
        'orjson>=3.0.2',
        'pyexcel>=0.6.4',
        'pandas>=1.0.3',
        'numpy>=1.18.1',
        'textdistance>=4.1.5',
        'databases>=0.3.2',
        # 'neo4j>=4.0.1',
        'rich>=9.5.0',
        'orm>=0.1.5',
        'scikit-learn>=0.22',
        'python-slugify>=4.0.0',
        'cachetools>=4.1.0',
        'click>=7.1.2'
     ],
    license="MIT",
    author_email="1730416009@stu.suda.edu.cn",
    maintainer="ZeFeng Zhu",
    maintainer_email="1730416009@stu.suda.edu.cn",
    description="Profiling Protein Structures from Protein Data Bank and integrate various resources.🏄‍♂️",
    long_description=readme,
    long_description_content_type="text/markdown",
    url="https://github.com/NatureGeorge/pdb-profiling",
    python_requires=">=3.6.*",
    classifiers=[
        "Operating System :: OS Independent",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9"
    ],
    )

"""
Packaged By
python setup.py sdist bdist_wheel
twine upload --repository-url https://test.pypi.org/legacy/ dist/*
twine upload --repository-url https://upload.pypi.org/legacy/ dist/*

> https://packaging.python.org/tutorials/packaging-projects/
"""
