Installation
------------

There are 4 methods of installation available (choose one):

Through the ``conda`` package manager (`Anaconda Cloud <https://anaconda.org/conda-forge/bamnostic>`__)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: bash

    # first, add the conda-forge channel to your conda build
    conda config --add channels conda-forge

    # now bamnostic is available for install
    conda install bamnostic

Through the Python Package Index (`PyPI <https://pypi.org/>`__)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: bash

    pip install bamnostic

    # or, if you don't have superuser access
    pip install --user bamnostic

Through pip+Github
~~~~~~~~~~~~~~~~~~

.. code:: bash

    # again, use --user if you don't have superuser access
    pip install -e git+https://github.com/betteridiot/bamnostic.git

    # or, if you don't have superuser access
    pip install --user -e git+https://github.com/betteridiot/bamnostic.git

Traditional GitHub clone
~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: bash

    git clone https://github.com/betteridiot/bamnostic.git
    cd bamnostic
    pip install -e .

    # or, if you don't have superuser access
    pip install --user -e .

