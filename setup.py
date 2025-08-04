from setuptools import setup, find_packages

setup(
    name='qm9nmr-ml',
    version='0.0.0',
    packages=find_packages(),
    package_data={'qm9nmr-ml': ['data/*']},
    author='Raghunathan Ramakrishnan',
    author_email='raghu.rama.chem@gmail.com',
    url='https://github.com/moldis-group/qm9nmr-ml',
    license='MIT License',
    description='qm9nmr-ml: A Python-based ML model trained on QM9NMR for 13C chemical shift predictions.',
    long_desc_type="text/markdown",
    install_requires=[ 'pandas', 'numpy','matplotlib' ],
    include_package_data=True
)

