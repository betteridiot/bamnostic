from setuptools import setup

setup(
    name='bamnostic',
    version='0.1',
    description='Pure Python, OS-agnostic Binary Alignment Map (BAM) random access and parsing tool',
    url='https://github.com/betteridiot/bamnostic/',
    author='Marcus D. Sherman',
    author_email='mdsherm@umich.edu',
    license='BSD 3-Clause',
    packages=['bamnostic', 'tests'],
    zip_safe=False
)
