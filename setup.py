try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

with open('README.md', 'r', encoding='utf-8') as fh:
    long_description = fh.read()

setup(
    name='l-giremi',
    version='0.1.8',
    author='Zhiheng Liu',
    author_email='wolfsonliu@live.com',
    description='a software for analysis of RNA editing sites from long-read RNA-seq data',
    long_description=long_description,
    long_description_content_type='text/markdown',
    url='https://github.com/xiaolab/L-GIREMI',
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)',
        'Operating System :: OS Independent',
        'Topic :: Scientific/Engineering :: Bio-Informatics'
    ],
    install_requires=[
        'pysam', 'numpy', 'pandas', 'sklearn'
    ],
    scripts = [
        'bin/calculate_site_splice_mi',
        'bin/correct_splice_site',
        'bin/get_aei',
        'bin/get_read_intron',
        'bin/get_read_mismatch',
        'bin/get_read_site',
        'bin/get_read_splice',
        'bin/split_bam_by_site',
        'l-giremi'
    ],
    python_requires='>=3.6',
)
