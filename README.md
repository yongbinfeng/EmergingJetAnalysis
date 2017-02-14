<div id="table-of-contents">
<h2>Table of Contents</h2>
<div id="text-table-of-contents">
<ul>
<li><a href="#orgf23ff69">1. EmergingJetAnalysis</a>
<ul>
<li><a href="#org524970b">1.1. Sub-packages</a></li>
<li><a href="#orgcd09575">1.2. Dependencies</a>
<ul>
<li><a href="#org07c449b">1.2.1. Cog - </a></li>
</ul>
</li>
<li><a href="#orgb1633ac">1.3. Quick instructions</a>
<ul>
<li><a href="#orgbb93360">1.3.1. To checkout latest code</a></li>
<li><a href="#org44bd7b4">1.3.2. To setup environment and compile</a></li>
<li><a href="#org1c6dcb9">1.3.3. To run main example config file</a></li>
</ul>
</li>
<li><a href="#orge11c53a">1.4. General steps</a></li>
<li><a href="#org7d19470">1.5. Example config file</a></li>
</ul>
</li>
</ul>
</div>
</div>

<a id="orgf23ff69"></a>

# EmergingJetAnalysis

Main repository for Emerging Jet analysis


<a id="org524970b"></a>

## Sub-packages

-   EmJetAnalyzer: Main analyzer - analyzes input jets and outputs flat NTuple
-   JetFilter: Applies basic kinematic selection for the signal region and outputs selected jets for analysis
-   WJetFilter: Applies W+jet selection and outputs selected jets for analysis


<a id="orgcd09575"></a>

## Dependencies


<a id="org07c449b"></a>

### Cog - <http://nedbatchelder.com/code/cog/index.html>

Cog is used for code generation of classes with repetitive structure, e.g. `EmergingJetAnalyzer/interface/OutputTree.h`
To install, simply execute:

    pip install --user cogapp


<a id="orgb1633ac"></a>

## Quick instructions


<a id="orgbb93360"></a>

### To checkout latest code

    cmsrel CMSSW_7_6_3
    cd CMSSW_7_6_3/src
    cmsenv
    git cms-merge-topic yhshin11:fix-GetTrackTrajInfo # Needed for track trajectory info.
    git clone https://github.com/yhshin11/EmergingJetAnalysis.git


<a id="org44bd7b4"></a>

### To setup environment and compile

    cmsenv
    cd EmergingJetAnalysis
    ./buildAll.sh


<a id="org1c6dcb9"></a>

### To run main example config file

    cmsRun Configuration/test/test_cfg.py


<a id="orge11c53a"></a>

## General steps

For our analysis involving Emerging Jets, we follow the following steps:

1.  Skim the data for useful events and select jets (JetFilter/WJetFilter)
2.  Run the analyzer over the selected events/jets and produce ROOT files containing simple TTree objects (EmJetAnalyzer)
3.  Make histograms/plots from the flat ROOT files


<a id="org7d19470"></a>

## Example config file

An example config file is provided at `Configuration/test/test_cfg.py`. Make a copy of this file and edit `process.source` to run over an appropriate EDM file (AOD/AODSIM). You can now run the config file using cmsRun.

The config file comes with some command line options (defined at the top of the file):

-   `data` (default: `0`) - Must be set to 1 when running over real data and not MC.
-   `steps` (default: `skim,analyze`) - Determines which of the provided steps to execute. Setting `skim` turns on the basic selection step. Setting `analyze` turns on the analyzer.
-   `sample` (default: `signal`) - This sets the module parameters to be appropriate for different samples, e.g. if `sample=wjet` is set, the W+jet kinematic selection (WJetFilter) is used instead of the signal region kinematic selection (JetFilter). The details of the configuration are defined in `Configuration/python/emJetTools.py`. (Warning: The logic here may change while the analysis is ongoing.)

For example to run both the skim and analyze steps over signal MC, simply do:

    cmsRun Configuration/test/test_cfg.py data=0 steps=skim,analyze sample=signal

Since the provided options are exactly the default, you could also simply do:

    cmsRun Configuration/test/test_cfg.py

Sometimes you may want to run the two steps separately. For instance since the skimming step is more stable while the analyze step is more in flux, you may want to run only the skimming step over a large dataset, so that you can run different versions of the analyze step later on. To run only the skimming step over data events, you can do:

    cmsRun Configuration/test/test_cfg.py data=1 steps=skim

Note: By default, EmJetAnalyzer expects as input, jets that were selected by WJetFilter or JetFilter. This means that you cannot only run the analyze step over a bare AOD/AODSIM file.

