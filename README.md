# Yade-interFoam coupling
This is a development initiated at OceanFOAM2024 hackathon.

In this repo we aim to couple Yade with interFoam solver. 

Here we follow the paradigm and implementations found [here](https://gitlab.com/puigmontella/yade-trunk/-/tree/newofoamVersions_checks/examples/openfoam/example_pimpleFoamYade?ref_type=heads) 

# Compile

Assuming you have downloaded and installed Yade you should 1) copy the folder FoamYade_two_phase to /pkg/openfoam/coupling/, 2) copy the folder interFoamYADEv2312 to /pkg/openfoam/coupling/Solvers, 3) compile the FoamYade_two_phase and interFoamYADEv2312 using wmake in the respective folders.

# Tutorial

There is a simple tutorial of a half filled unit cube with one light particle incerded under the air-water interface

# Status 

This is work in progress....