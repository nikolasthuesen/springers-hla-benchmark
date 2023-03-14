import pathlib
from setuptools import setup, find_packages

file_folder = pathlib.Path(__file__).parent

# get the dependencies and installs
with (file_folder / 'requirements.txt').open() as f:
    all_reqs = f.read().split('\n')

setup(
    name='hla_typing_benchmark',
    version='0',
    license='Proprietary',
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Developers',
        'Programming Language :: Python :: 3',
    ],
    keywords='',
    packages=find_packages(where='src', exclude=['docs', 'tests*']),
    package_dir={'':'src'},
    include_package_data=True,
    author='Nikolas Thuesen',
    install_requires=all_reqs,
    setup_requires=['wheel'],
    entry_points={
            'console_scripts': [
                'parse_typing_results = hla_typing_benchmark.parse_typing_results:main',
                'create_gold_standard = hla_typing_benchmark.create_gold_standard:main',
                'summarise_results = hla_typing_benchmark.summarise_results:main',
            ]
        },
)
