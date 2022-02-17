# bioinformatics_project_2_2020

Repo for the second project of the Bioinformatics & Network Medicine 2020 at Sapienza. To install all dependencies run
```shell
$ pip3 -r requirements.txt
```

## Notes

### Fixing broken EDF headers

Before using EDF files fix the header using EDFBrowser 
```shell
$ sudo apt install edfbrowser
```
See differences in the screenshot.

### Installing bctpy

To install bctpy for python

```shell  
$ git clone https://github.com/aestrivex/bctpy
$ cd bctpy
$ python3 setup.py build
$ python3 setup.py install
$ sudo mv motif34lib.mat /usr/local/lib/python3.6/dist-packages/bctpy-0.5.2-py3.6.egg/bct
$ cd \usr\local\lib\python3.6\dist-packages\bctpy-0.5.2-py3.6.egg\bct\algorithms
```

Then open file `motifs.py` and change array sizes to `200` at lines 701-702 (refer to the second screenshot). 

## Authors
- [Giammarco D'Alessandro](https://github.com/giamdalessandro)
- [Luca Gioffrè](https://github.com/lukfre)
- Narges Sharghi
