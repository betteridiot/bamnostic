|Documentation Status| |Conda Version| |PyPI version| |Maintainability|

|status| |DOI| |License|

+----------+-------------------------+
| Platform | Build Status            |
+==========+=========================+
| Linux    | |Build Status TravisCI| |
+----------+-------------------------+
| Windows  | |Build status Appveyor| |
+----------+-------------------------+
| conda    | |noarch|                |
+----------+-------------------------+

+-------+-------------------+
| Host  | Downloads         |
+=======+===================+
| PyPI  | |Downloads|       |
+-------+-------------------+
| conda | |Conda Downloads| |
+-------+-------------------+

BAMnostic
=========

a *pure Python*, **OS-agnositic** Binary Alignment Map (BAM) file parser
and random access tool.

Note:
~~~~~

Documentation can be found at `here`_ or by going to this address:
http://bamnostic.readthedocs.io. Documentation was made available
through `Read the Docs`_.

--------------

Installation
------------

There are 4 methods of installation available (choose one):

Through the ``conda`` package manager (`Anaconda Cloud`_)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: bash

   # first, add the conda-forge channel to your conda build
   conda config --add channels conda-forge

   # now bamnostic is available for install
   conda install bamnostic

Through the Python Package Index (`PyPI`_)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: bash

   pip install bamnostic

   # or, if you don't have superuser access
   pip install --user bamnostic

Through pip+Github
~~~~~~~~~~~~~~~~~~

.. code:: bash

   # again, use --user if you don't have superuser access
   pip install -e git+https://github.com/betteridiot/bamnostic.git#egg=bamnostic

   # or, if you don't have superuser access
   pip install --user -e git+https://github.com/betteridiot/bamnostic.git#bamnostic#egg=bamnostic

Traditiona
~~~~~~~~~~

.. _here: http://bamnostic.readthedocs.io/en/latest/
.. _Read the Docs: https://readthedocs.org/
.. _Anaconda Cloud: https://anaconda.org/conda-forge/bamnostic
.. _PyPI: https://pypi.org/

.. |Documentation Status| image:: https://readthedocs.org/projects/bamnostic/badge/?version=latest
   :target: https://bamnostic.readthedocs.io/en/latest/?badge=latest
.. |Conda Version| image:: https://img.shields.io/conda/vn/conda-forge/bamnostic.svg
   :target: https://anaconda.org/conda-forge/bamnostic
.. |PyPI version| image:: https://badge.fury.io/py/bamnostic.svg
   :target: https://badge.fury.io/py/bamnostic
.. |Maintainability| image:: https://api.codeclimate.com/v1/badges/d7e36e72f109c598c86d/maintainability
   :target: https://codeclimate.com/github/betteridiot/bamnostic/maintainability
.. |status| image:: http://joss.theoj.org/papers/9952b35bbb30ca6c01e6a27b80006bd8/status.svg
   :target: http://joss.theoj.org/papers/9952b35bbb30ca6c01e6a27b80006bd8
.. |DOI| image:: https://zenodo.org/badge/DOI/10.5281/zenodo.1341959.svg
   :target: https://doi.org/10.5281/zenodo.1341959
.. |License| image:: https://img.shields.io/badge/License-BSD%203--Clause-blue.svg
   :target: https://github.com/betteridiot/bamnostic/blob/master/LICENSE
.. |Build Status TravisCI| image:: https://travis-ci.org/betteridiot/bamnostic.svg?branch=master
   :target: https://travis-ci.org/betteridiot/bamnostic
.. |Build status Appveyor| image:: https://ci.appveyor.com/api/projects/status/y95q02gkv3lgmlf4/branch/master?svg=true
   :target: https://ci.appveyor.com/project/betteridiot/bamnostic/branch/master
.. |noarch| image:: https://img.shields.io/circleci/project/github/conda-forge/bamnostic-feedstock/master.svg?label=noarch
   :target: https://circleci.com/gh/conda-forge/bamnostic-feedstock
.. |Downloads| image:: http://pepy.tech/badge/bamnostic
   :target: http://pepy.tech/project/bamnostic
.. |Conda Downloads| image:: https://img.shields.io/conda/dn/conda-forge/bamnostic.svg
   :target: https://anaconda.org/conda-forge/bamnostic