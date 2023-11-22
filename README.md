# Modules for `nanoAOD-tools`

This repo provides some modules that can be used with the
[`nanoAOD-tools`](https://github.com/cms-nanoAOD/nanoAOD-tools)
framework for post-processing nanoAOD.

The modules add new branches with like lepton scale factors, using the
[`correctionlib`](https://github.com/cms-nanoAOD/correctionlib)
tool, and official
[JSON correction files](https://gitlab.cern.ch/cms-nanoAOD/jsonpog-integration)
provided by the CMS POGs.


## Installation

### Install CMSSW
First, install your favorite CMSSW release.
If you rely on `correctionlib` and `python3`,
it is strongly recommend to use CMSSW 11.3 or newer.
For example,
<table>
<tr>
<td>
CMSSW 11.3.4 ([Combine v9](http://cms-analysis.github.io/HiggsAnalysis-CombinedLimit/))
</td>
<td>
CMSSW 12.4.8
</td>
</tr>
<tr>
<td>

```bash
export SCRAM_ARCH=slc7_amd64_gcc900
cmsrel CMSSW_11_3_4
cd CMSSW_11_3_4/src
cmsenv
```
</td>
<td>

```bash
export SCRAM_ARCH=slc7_amd64_gcc10
cmsrel CMSSW_12_4_8
cd CMSSW_12_4_8/src/
cmsenv
```
</td>
</tr>
</table>



### Install `nanoAOD-tools` (for CMSSW 12.6 or older)
Install [`nanoAOD-tools`](https://github.com/cms-nanoAOD/nanoAOD-tools) to process nanoAOD files.
Note that starting from CMSSW 13.0, a basic version of `nanoAOD-tools` is included.
To install the standalone version, please do
```bash
cd $CMSSW_BASE/src/
git clone https://github.com/cms-nanoAOD/nanoAOD-tools.git PhysicsTools/NanoAODTools
```

### Install `NATModules`
Now install this repository as `PhysicsTools/NATModules` and compile (build) everything:
```bash
cd $CMSSW_BASE/src/
cmsenv
git clone git@github.com:cms-cat/nanoAOD-tools-modules.git PhysicsTools/NATModules
scram b
```

### Install `correctionlib` (for CMSSW 11.2 or older)
Starting from CMSSW 11.3,
[`correctionlib`](https://github.com/cms-nanoAOD/correctionlib)
should come pre-installed.
To install yourself for older CMSSW versions,
please see have a look at [documentation](https://cms-nanoaod.github.io/correctionlib/).
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
Alternatively, this repository is regularly synchronized to `/cvmfs/`,
so if your system has access, you can copy the latest version
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
