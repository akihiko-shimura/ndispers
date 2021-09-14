from setuptools import setup, find_packages
from codecs import open
from os import path

root_dir = path.abspath(path.dirname(__file__))

def _requirements():
    return [name.rstrip() for name in open(path.join(root_dir, 'requirements.txt')).readlines()]

with open('README.md', encoding='utf-8') as f:
    long_description = f.read()

setup(
    name="ndispers",
    version='0.1.8',
    packages=find_packages(),

    license='MIT License',
    author='Akihiko Shimura',
    author_email='shimura.akihiko@opt-oxide.com',
    url='https://github.com/akihiko-shimura/ndispers',
    description='A Python package for calculating refractive index dispersion of various materials.',
    long_description=long_description,

    install_requires=_requirements(),

    keywords='nonlinear optics physics',
    classifiers=[
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
    ],
)
