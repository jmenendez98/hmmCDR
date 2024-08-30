from setuptools import setup, find_packages

setup(
    name='hmmCDR',
    version='0.1.0',
    description="Find CDR locations using bedmethyl file and CenSat annotations.",
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    url="https://github.com/jmenendez98/hmmCDR",
    author='Julian Menendez',
    author_email='jmmenend@ucsc.edu',
    license="MIT",
    packages=find_packages(),
    install_requires=[
        'numpy>=1.21.5',
        'pandas>=1.3.5',
        'pybedtools>=0.8.1',
        'hmmlearn>=0.3.0'
    ],
    entry_points={
        'console_scripts': [
            'hmmCDR=hmmCDR.hmmCDR:main',
            'hmmCDRprior=hmmCDR.hmmCDRprior:main'
        ]
    },
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
    ],
    python_requires='>=3.7',
)