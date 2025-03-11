from setuptools import setup, find_packages

setup(
    name='multiplexdesigner',
    version='0.1.0',
    packages=find_packages(),
    install_requires=[
        # List your package dependencies here
    ],
    author='Stefan Filges Name',
    author_email='stefan.filges@pm.me',
    description='A package for multiplex PCR assay design',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    url='https://github.com/sfilges/multiplex_designer',
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
    ],
    python_requires='>=3.6',
)