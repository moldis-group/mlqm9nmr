from setuptools import setup, find_packages

setup(
    name='mlqm9nmr',
    version='0.0.0',
    packages=find_packages(),
    package_data={'mlqm9nmr': ['data/*']},
    author='Raghunathan Ramakrishnan',
    author_email='raghu.rama.chem@gmail.com',
    url='https://github.com/moldis-group/mlqm9nmr',
    license='MIT License',
    description='mlqm9nmr: A Python-based ML model trained on QM9NMR for 13C chemical shift predictions.',
    long_desc_type="text/markdown",
    install_requires=[ 'pandas', 'numpy','matplotlib' ],
    include_package_data=True
)

