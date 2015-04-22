# xyz2top Python Module

<!---
Wercker badge
-->

This module provides generates genuine topology files from a molecular geometry provided in [XYZ format](http://en.wikipedia.org/wiki/XYZ_file_format). Building such a file requires some specific and arbitrary definition of what a molecular bond is and how it is detected.

We will try to summarize here the decisions that leads to our final topology file:

- This module only attemps to identify covalent bonds
- covalent bonds are detected as such when the distance between 2 atoms is simply smaller or equal to the sum of their respective [covalent radius](http://chemwiki.ucdavis.edu/Theoretical_Chemistry/Chemical_Bonding/General_Principles/Covalent_Bond_Distance,_Radius_and_van_der_Waals_Radius).

  ![test text](http://chemwiki.ucdavis.edu/@api/deki/files/10118/covalent_vanderwaals.png)

http://chemwiki.ucdavis.edu/Theoretical_Chemistry/Chemical_Bonding/General_Principles/Covalent_Bond_Distance,_Radius_and_van_der_Waals_Radius
The sum of the covalent radii can be optionnally multiplied by a coefficient (by default equal to 1.3. With that value, our module reproduces the bonds detected in a few tested molecules in Avogadro)) before being compared to the actual distance of the pair of atoms (e.g. if the coefficient is increased, potentially more pairs of atoms will be detected as being connected by a covalent bond).


## Table of contents

- [Quick start](#quick-start)
- [Testing](#testing)
- [Dependencies](#dependencies)
- [License](#license)
- [Issues](#issues)


## Quick start

You can install the xyz2top module simply with the following commands: 
```shell
git clone https://github.com/Patechoc/xyz2top.git
cd xyz2top
sudo python setup.py install
```

You can install the module locally with a virtual environment and the dependencies that were tested last for passing the tests:
```shell
git clone https://github.com/Patechoc/xyz2top.git
cd xyz2top
sudo easy_install pip
virtualenv myPackage
source myPackage/bin/activate
sudo python setup.py install
py.test tests/
deactivate
```

or installing the development version:

```shell
git clone https://github.com/Patechoc/xyz2top.git
cd xyz2top
sudo easy_install pip
virtualenv myPackage
source myPackage/bin/activate
sudo ./myPackage/bin/pip install -r requirements-dev.txt
py.test tests/
deactivate
```



## Testing

```shell
cd tests
py.test 
```

or some specific test module:

```shell
cd tests
py.test -q test_main.py  # to run a single test file
```



## Dependencies

- elements.py is an external module containing information of the periodic table, among which some experimental covalent radii to atoms used for identifying covalent bond according to our simple model of such bonds.

## License

The code which makes up this Python project template is licensed under the MIT/X11 license. Feel free to use it in your free software/open-source or proprietary projects.


## Issues

Please report any bugs or requests that you have using the GitHub issue tracker!
