# Channel-noMOC-Darwin

[![DOI](https://zenodo.org/badge/194902741.svg)](https://zenodo.org/badge/latestdoi/194902741)

Repository for the Darwin package and configuration files for the MITgcm channel model with surface seasonal forcings but without a northern sponge layer, i.e. no overturning circulation.
The [Darwin package used](/pkg/darwin) is not compatible with the most up-to-date version of MITgcm. To revert to the specific version, run:

`git clone https://github.com/MITgcm/MITgcm.git`
`cd MITgcm`
`git checkout checkpoint63n.`

The config files to run the Darwin package are [here](/Darwin-config/code) and namelist files [here](/Darwin-config/run_proper_phase/run).

Jupyter notebooks to generate the input files are [here](/Darwin-config/notebooks).
