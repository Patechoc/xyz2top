=========================
 xyz2top Python Module
=========================

<!---
Wercker badge
-->

This module provides builds some kind of topology from a molecular geometry provided in [XYZ format](http://en.wikipedia.org/wiki/XYZ_file_format). Building such a file requires some specific and arbitrary definition of what a molecular bond is and how it is detected.

We will try to summarize here the decisions that leads to our final topology file:
- criteria #1
- criteria #2
- criteria #3



## Table of contents

- [Quick start](#quick-start)
- [Testing](#testing)
- [License](#license)
- [Issues](#issues)

## Quick start

You can install the xyz2top module simply with the following commands: 
```
git clone https://github.com/Patechoc/xyz2top.git
cd xyz2top
sudo python setup.py install
```

You can install the module locally with a virtual environment and the dependencies that were tested last for passing the tests:
```
git clone https://github.com/Patechoc/xyz2top.git
cd xyz2top
sudo easy_install pip
virtualenv myPackage
source myPackage/bin/activate
sudo ./myPackage/bin/pip install -r requirements.txt
cd test; python -m unittest discover -v  ## optional but quick ;)
deactivate
```


## Testing

```
cd test
python -m unittest discover -v
```


## License

The code which makes up this Python project template is licensed under the MIT/X11 license. Feel free to use it in your free software/open-source or proprietary projects.


## Issues

Please report any bugs or requests that you have using the GitHub issue tracker!
