from setuptools import setup, find_packages

def main():
    setup(
        name='ndispers',
        version='0.1.3',
        packages=find_packages(exclude=['tests']),

        author='Akihiko Shimura',
        author_email='shimura.akihiko@opt-oxide.com',

        url='https://github.com/akihiko-shimura/ndispers',
        description='A Python package for calculating refractive index dispersion of various materials.',
        long_description=open('README.md').read(),
        long_description_content_type='text/markdown',

        python_requires='~=3',

        classifiers=[
            'License :: OSI Approved :: MIT License',
            'Programming Language :: Python :: 3',
            'Programming Language :: Python :: 3.9',
        ],
        )

if __name__ == '__main__':
    main()
