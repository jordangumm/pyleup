from setuptools import setup, find_packages

setup(
    name='pyleup',
    version='0.1',
    packages=find_packages(),
    scripts=['bin/pyleup'],

    install_requires=[],

    author='Jordan Gumm',
    author_email='jordan@variantanalytics.com',
    description='A tool for creating visualizations for NGS reference read pileups'
)
