from setuptools import setup
import os
import sys


def readme():
    path = os.path.dirname(__file__) or "."
    with open(os.path.join(path, "README.rst")) as rst:
        return rst.read()


setup(
    name="bamnostic",
    version=open("version").read().strip(),
    description="Pure Python, OS-agnostic Binary Alignment Map (BAM) random access and parsing tool",
    long_description=readme(),
    url="https://github.com/betteridiot/bamnostic/",
    author="Marcus D. Sherman",
    author_email="mdsherman@betteridiot.tech",
    license="BSD 3-Clause",
    extras_require={'test': ['pytest']},
    packages=["bamnostic", "tests"],
    package_dir={"bamnostic": "./bamnostic", "tests": "./tests"},
    package_data={
        "bamnostic": [
            "data/*",
            "LICENSE",
            "CONTRIBUTING.md",
            "CODE_OF_CONDUCT.md",
            "version",
        ]
    },
    classifiers=[
        "Development Status :: 5 - Production/Stable",
        "Intended Audience :: Developers",
        "Intended Audience :: End Users/Desktop",
        "Intended Audience :: Science/Research",
        "Intended Audience :: Information Technology",
        "License :: OSI Approved :: BSD License",
        "Natural Language :: English",
        "Operating System :: Unix",
        "Operating System :: Microsoft :: Windows",
        "Operating System :: MacOS",
        "Programming Language :: Python :: 2.7",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.0",
        "Programming Language :: Python :: 3.1",
        "Programming Language :: Python :: 3.2",
        "Programming Language :: Python :: 3.3",
        "Programming Language :: Python :: 3.4",
        "Programming Language :: Python :: 3.5",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Programming Language :: Python :: 3.12",
        "Programming Language :: Python :: 3.13",
        "Programming Language :: Python :: 3.14",
        "Programming Language :: Python :: 3.15",
        "Programming Language :: Python :: Implementation :: CPython",
        "Programming Language :: Python :: Implementation :: PyPy",
        "Topic :: Scientific/Engineering",
        "Topic :: Scientific/Engineering :: Information Analysis",
    ],
    keywords="BAM pysam genomics genetics struct",
    include_package_data=True,
    zip_safe=False,
)
