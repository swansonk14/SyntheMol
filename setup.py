from pathlib import Path
from setuptools import find_packages, setup

# Load version number
__version__ = ''
version_file = Path(__file__).parent.absolute() / 'synthemol' / '_version.py'

with open(version_file) as fd:
    exec(fd.read())

# Load README
with open('README.md', encoding='utf-8') as f:
    long_description = f.read()

setup(
    name='synthemol',
    version=__version__,
    author='Kyle Swanson',
    author_email='swansonk.14@gmail.com',
    description='synthemol',
    long_description=long_description,
    long_description_content_type='text/markdown',
    url='https://github.com/swansonk14/SyntheMol',
    download_url=f'https://github.com/swansonk14/SyntheMol/archive/refs/tags/v_{__version__}.tar.gz',
    license='MIT',
    packages=find_packages(),
    package_data={'synthemol': ['py.typed', 'resources/**/*']},
    entry_points={
        'console_scripts': [
            'synthemol=synthemol.generate.generate:generate_command_line'
        ]
    },
    install_requires=[
        'chemfunc',
        'chemprop',
        'descriptastorus',
        'matplotlib',
        'numpy',
        'pandas',
        'torch',
        'rdkit',
        'scikit-learn',
        'scipy',
        'tqdm',
        'typed-argument-parser>=1.8.0'
    ],
    python_requires='>=3.10',
    classifiers=[
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.10',
        'Programming Language :: Python :: 3.11',
        'Programming Language :: Python :: 3.12',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
        'Typing :: Typed'
    ],
    keywords=[
        'machine learning',
        'drug design',
        'generative models',
        'synthesizable molecules'
    ]
)
