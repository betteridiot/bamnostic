from setuptools.command.test import test as TestCommand
from setuptools import setup
import os
import sys


class PyTest(TestCommand):
    user_args = [('pytest-args=', 'a', 'Arguments to pass to py.test')]

    def initialize_options(self):
        TestCommand.initialize_options(self)
        self.pytest_args = []

    def run_tests(self):
        import pytest
        errno = pytest.main(self.pytest_args)
        if errno:
            sys.exit(errno)
        else:
            errno = pytest.main(['--doctest-modules', '-Wignore'])
            sys.exit(errno)


def readme():
    path = os.path.dirname(__file__) if os.path.dirname(__file__) else '.'
    with open(path + '/README.rst') as rst:
        return rst.read()


setup(
    name='bamnostic',
    version=open('version').read(),
    description='Pure Python, OS-agnostic Binary Alignment Map (BAM) random access and parsing tool',
    long_description=readme(),
    url='https://github.com/betteridiot/bamnostic/',
    author='Marcus D. Sherman',
    author_email='mdsherm@umich.edu',
    license='BSD 3-Clause',
    # setup_requires=['pytest-runner'],
    tests_require=['pytest'],
    cmdclass = {'test' : PyTest},
    packages=['bamnostic', 'tests'],
    package_dir={'bamnostic': './bamnostic', 'tests': './tests'},
    package_data={'bamnostic': ['data/*', 'LICENSE', 'CONTRIBUTING.md', 'CODE_OF_CONDUCT.md', 'version']},
    classifiers=[
        'Development Status :: 5 - Production/Stable',
        'Intended Audience :: Developers',
        'Intended Audience :: End Users/Desktop',
        'Intended Audience :: Science/Research',
        'Intended Audience :: Information Technology',
        'License :: OSI Approved :: BSD License',
        'Natural Language :: English',
        'Operating System :: Unix',
        'Operating System :: Microsoft :: Windows',
        'Operating System :: MacOS',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.0',
        'Programming Language :: Python :: 3.1',
        'Programming Language :: Python :: 3.2',
        'Programming Language :: Python :: 3.3',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: Implementation :: CPython',
        'Programming Language :: Python :: Implementation :: PyPy',
        'Topic :: Scientific/Engineering',
        'Topic :: Scientific/Engineering :: Information Analysis'
    ],
    keywords='BAM pysam genomics genetics struct',
    include_package_data=True,
    zip_safe=False
)
