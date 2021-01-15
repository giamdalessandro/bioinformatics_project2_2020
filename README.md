# bioinformatics_project_2_2020

Repo for the second project of the Bioinformatics & Network Medicine 2020 at Sapienza

## notes

Before using EDF files fix the header using EDFBrowser `$ sudo apt install edfbrowser`.
See differences in the screenshot.

### bctpy
To install bctpy for python

```shell  
$ git clone https://github.com/aestrivex/bctpy
$ cd bctpy
$ python3 setup.py build
$ python3 setup.py install
$ sudo mv motif34lib.mat /usr/local/lib/python3.6/dist-packages/bctpy-0.5.2-py3.6.egg/bct
```