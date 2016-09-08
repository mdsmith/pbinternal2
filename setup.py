import os
from setuptools import setup, find_packages


version = __import__('pbinternal2').get_version()

_REQUIREMENTS_FILE = 'REQUIREMENTS.txt'


def _get_local_file(file_name):
    return os.path.join(os.path.dirname(__file__), file_name)


def _get_requirements(file_name):
    with open(file_name, 'r') as f:
        lines = f.readlines()
    reqs = [l for l in lines if not l.startswith("#")]
    return reqs

setup(
    name='pbinternal2',
    version=version,
    author='pbiDevNet',
    author_email='pbiDevNet@pacificbiosciences.com',
    license='LICENSE.txt',
    packages=find_packages(),
    zip_safe=False,
    install_requires=_get_requirements(_get_local_file(_REQUIREMENTS_FILE)),
    entry_points={'console_scripts': [
        'makeMissingAdapterReport = pbinternal2.report.missing_adapter:main',
        'makeRainbowReport = pbinternal2.report.rainbow:main',
        'findAdapterIssues = pbinternal2.util.find_adapter_issues:main'
        ]}
)
