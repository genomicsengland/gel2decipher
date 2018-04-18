from distutils.core import setup
from setuptools import find_packages

setup(
    name='gel2decipher',
    version='0.1.0',
    packages=find_packages(),
    scripts=['scripts/gel2decipher.py'],
    url='',
    license='',
    author='priesgo',
    author_email='pablo.ferreiro@genomicsengland.co.uk',
    description='',
    install_requires=[
        'requests',
        'pyyaml',
        'gelreportmodels==6.0.6',
        'booby==0.7.0'
    ]
)