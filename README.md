# Modules for `NanoAODTools`

This repo provides some modules that can be used with the
[`NanoAODTools`](https://github.com/cms-sw/cmssw/tree/master/PhysicsTools/NanoAODTools)
framework for post-processing CMS's nanoAOD files.

The modules add new branches like lepton scale factors, using the
[`correctionlib`](https://github.com/cms-nanoAOD/correctionlib)
tool, and official
[JSON correction files](https://gitlab.cern.ch/cms-nanoAOD/jsonpog-integration)
provided by the CMS POGs.


## Installation

### Install CMSSW
First, install your favorite CMSSW release. Recent correction files require correctionlib >=2.6 which is provided in recent CMSSW releases, eg. CMSSW 14.

```bash
cmsrel CMSSW_14_1_6
cd CMSSW_14_1_6/src
cmsenv
```

### Install `NATModules`
Now install this repository as `PhysicsTools/NATModules` and compile (build) everything:
```bash
cd $CMSSW_BASE/src/
cmsenv
git clone git@github.com:cms-cat/nanoAOD-tools-modules.git PhysicsTools/NATModules
scram b
```

### For very old CMSSW releases

In CMSSW 12 or older, 13_0_X with X < 16, 13_1_X with X < 5, 13_2_X, `NanoAODTools` must be installed manually (it is included in release starting from CMSSW 13_0_16, 13_1_5 and 13_3_0). It can be [taken from the CMSSW release area](CMSSW https://github.com/cms-sw/cmssw/tree/master/PhysicsTools/NanoAODTools).

In CMSSW 11.2 or older, [`correctionlib`](https://github.com/cms-nanoAOD/correctionlib) must also be installed manually (it is distrubuted with CMSSW starting from CMSSW 11.3)

To install it yourself, please see have a look at [documentation](https://cms-nanoaod.github.io/correctionlib/).
Note that `correctionlib` works best for `python3`.

### Install `correctionlib` data
Finally, install the `correctionlib` JSON files provided by CMS
into `PhysicsTools/NATModules/data` from the
[`cms-nanoAOD/jsonpog-integration` repository](https://gitlab.cern.ch/cms-nanoAOD/jsonpog-integration)
on GitLab:
```bash
cd $CMSSW_BASE/src/PhysicsTools/NATModules
git clone ssh://git@gitlab.cern.ch:7999/cms-nanoAOD/jsonpog-integration.git data
```
or clone via Kerberos, where `$USER` is your CERN lxplus name:
```bash
kinit $USER@CERN.CH
git clone https://$USER:@gitlab.cern.ch:8443/cms-nanoAOD/jsonpog-integration.git
```
Alternatively, this repository is regularly synchronized to `/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/`,
so if your system has access, you can copy (or link) the latest version
```bash
cd $CMSSW_BASE/src/PhysicsTools/NATModules
cp -r /cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration data
```


### Test the installation
Run a module, for example
```bash
cd $CMSSW_BASE/src/PhysicsTools/NATModules
python3 ./test/example_muonSF.py -i root://cms-xrd-global.cern.ch//store/mc/RunIISummer20UL16NanoAODv9/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/NANOAODSIM/20UL16JMENano_106X_mcRun2_asymptotic_v17-v1/2820000/11061525-9BB6-F441-9C12-4489135219B7.root
```


## Usage

To use in your own analysis, you can use the standalone scripts in [`test`](test/) as an example.
If you compiled this package correctly, you can import the modules in
[`PhysicsTools/NATModules/python/modules`](PhysicsTools/NATModules/python/modules)
as
```python
from PhysicsTools.NATModules.modules.muonSF import *
```
