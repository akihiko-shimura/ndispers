from setuptools import setup
from setuptools import find_packages

def main():
    setup(
        name='ndispers',
        version='0.1.2',
        zip_safe=False,
        packages=find_packages(),
        
        author='Akihiko Shimura'
        )

if __name__ == '__main__':
    main()
