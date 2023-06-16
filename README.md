# Brain network study during resting states

Second project of the Bioinformatics & Network Medicine 20/21 course, Sapienza University of Rome. 

We conducted the study of the a brain network during resting phase with eyes open and eyes closed. We used the EEG data of subject S003 from the dataset of PhysioNet (reachable [here](https://physionet.org/content/eegmmidb/1.0.0/)) to perform functional connectivity analysis. The study is divided into four steps; connectivity estimation and network creation, network analysis through local and global indices, motif analysis and community detection. We developed these tests in a python environment. Full report is avilable [here](./docs/report_proj_group2.pdf).


## Notes

### Run

To install all dependencies run
```shell
$ pip3 -r requirements.txt
```

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
- [Narges Sharghi](https://github.com/nargessharghi)
