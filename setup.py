from setuptools import setup, find_packages
from pathlib import Path
this_directory = Path(__file__).parent
long_description = (this_directory / "README.md").read_text()


setup(
    name='accfg',
    version='0.0.7',
    description='AccFG: Molecule functional group extraction and comparison',
    long_description=long_description,
    long_description_content_type='text/markdown',
    author='Xuan Liu',
    author_email='xliu254@illinois.edu',
    url='https://github.com/xuanliugit/AccFG',
    install_requires=[
        'requests',
        'importlib-metadata; python_version>"3.10"',
    ],
    packages=find_packages(
        include=['accfg*'], 
    ),
    include_package_data=True,
)