import pathlib
from setuptools import setup, find_packages

file_folder = pathlib.Path(__file__).parent

# get the dependencies and installs
with (file_folder / 'requirements.txt').open() as f:
    all_reqs = f.read().split('\n')

install_requires = [x.strip() for x in all_reqs if 'git+' not in x]
dependency_links = [x.strip().replace('git+', '') for x in all_reqs if x.startswith('git+')] #yapf: disable

setup(
    name='munin_gen0_designs',
    version='0',
    license='Proprietary',
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Developers',
        'Programming Language :: Python :: 3',
    ],
    keywords='',
    packages=find_packages(where='lib', exclude=['docs', 'tests*']),
    package_dir={'': 'lib'},
    include_package_data=True,
    author='Evaxion',
    install_requires=install_requires,
    dependency_links=dependency_links,
    setup_requires=['wheel'],
)
