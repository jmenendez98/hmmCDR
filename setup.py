# WIP: NEEDS UPDATING

from setuptools import setup, find_packages

setup(
    name='my_project',
    version='0.1.0',
    packages=find_packages(),
    install_requires=[
        'numpy',
        'pandas',
        'pybedtools'
    ],
    entry_points={
        'console_scripts': [
            'script1=my_package.script1:main',
            'script2=my_package.script2:main'
        ]
    },
    author='Julian Menendez',
    author_email='jmmenend@ucsc.edu',
    description='Uses a Hidden Markov Model to predict CDR sites',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    url='https://github.com/jmenendez98/hmmCDR',
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
    ],
    python_requires='>=3.6',
)