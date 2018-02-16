from setuptools import setup

setup(
    name='bamnostic',
    version='0.1',
    description='Pure Python, OS-agnostic Binary Alignment Map (BAM) indexing and parsing tool',
    url='https://github.com/betteridiot/bamnostic/',
    author='Marcus D. Sherman',
    author_email='mdsherm@umich.edu',
    license='CC BY 4.0',
    packages=['bamnostic','tests'],
    zip_safe=False)