xyz2top Python Module
=====================

.. raw:: html

   <!---
   Wercker badge
   -->

This module provides generates genuine topology files from a molecular
geometry provided in `XYZ
format <http://en.wikipedia.org/wiki/XYZ_file_format>`__. Building such
a file requires some specific and arbitrary definition of what a
molecular bond is and how it is detected.

We will try to summarize here the decisions that leads to our final
topology file:

-  This module only attemps to identify covalent bonds
-  covalent bonds are detected as such as when the distance between 2
   atoms is simply smaller or equal to the sum of their respective
   `covalent
   radius <http://chemwiki.ucdavis.edu/Theoretical_Chemistry/Chemical_Bonding/General_Principles/Covalent_Bond_Distance,_Radius_and_van_der_Waals_Radius>`__.
   The sum of the covalent radii can be optionnally multiplied by a
   coefficient (by default equal to 1) before being compared to the
   actual distance of the pair of atoms (e.g. if the coefficient is
   greater than 1., potentially more pairs of atoms will be detected as
   being connected by a covalent bond).

Table of contents
-----------------

-  `Quick start <#quick-start>`__
-  `Testing <#testing>`__
-  `Dependencies <#dependencies>`__
-  `License <#license>`__
-  `Issues <#issues>`__

Quick start
-----------

You can install the xyz2top module simply with the following commands:

::

    git clone https://github.com/Patechoc/xyz2top.git
    cd xyz2top
    sudo python setup.py install

You can install the module locally with a virtual environment and the
dependencies that were tested last for passing the tests:

::

    git clone https://github.com/Patechoc/xyz2top.git
    cd xyz2top
    sudo easy_install pip
    virtualenv myPackage
    source myPackage/bin/activate
    sudo ./myPackage/bin/pip install -r requirements-dev.txt
    cd tests; py.test
    deactivate

Testing
-------

::

    cd tests
    py.test -q test_main.py  # to run a single test file

or simply the following command for the full test suite:

::

    cd tests
    py.test 

Dependencies
------------

-  elements.py is an external module containing information of the
   periodic table, among which some experimental covalent radii to atoms
   used for identifying covalent bond according to our simple model of
   such bonds.

License
-------

The code which makes up this Python project template is licensed under
the MIT/X11 license. Feel free to use it in your free
software/open-source or proprietary projects.

Issues
------

Please report any bugs or requests that you have using the GitHub issue
tracker!
