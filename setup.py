from setuptools import setup, find_packages

setup(
    name='hmmCDR',
    version='0.3.0',
    description="Used pybedtools and hmmlearn to find precise locations of CDRs and subCDRs.",
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    url="https://github.com/jmenendez98/hmmCDR",
    author='Julian Menendez',
    author_email='jmmenend@ucsc.edu',
    license="MIT",
    packages=find_packages(),
    install_requires=[
        'numpy',
        'pandas',
        'pybedtools'
    ],
    entry_points={
        'console_scripts': [
            'hmmCDR=hmmCDR.hmmCDR:main',
            'hmmCDRprior=hmmCDR.hmmCDRprior:main',
            'hmmCDRparse=hmmCDR.hmmCDRparse:main'
        ]
    },
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
    ],
    python_requires='>=3.6',
)