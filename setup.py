from setuptools import setup, find_packages
from os import path

here = path.abspath(path.dirname(__file__))

# Get the long description from the README file
with open(path.join(here, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

setup(
    name="plex_stats",
    version="1.0.0",
    description='Calculates multiplex primer stats.',
    install_requires=[
        'Click',
    ],
    packages=['primer_stats'],
    entry_points='''
        [console_scripts]
        plex_stats=plex_stats.__main__:cli
    ''',
    long_description=long_description,
    long_description_content_type='text/markdown',
    url='https://github.com/tfursten/plex_stats',
    author='Tara Furstenau',
    author_email='tara.furstenau@nau.edu',
    classifiers=[ 
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
    ], 
    python_requires='>=3.6'
)
