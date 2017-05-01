# PyQMC

Version 0.01

Date: 2017-04-27

This is a Python computational toolbox for ab initio
quantum chemistry and condensed matter physics research.
It was originally written to support auxiliary-field quantum Monte
Carlo (AFQMC) calculations.
For AFQMC method, please refer to papers cited below.

This toolbox contains varied tools:

* output parsers for abinit, nwchem, gamess

* Gaussian basis tools & library

* matrix generator/handling for Hamiltonian & Slater determinants

* PWQMC / GAFQMC measurement tools (conversion to HDF5,
  reblocking tools, ...)

* site-specific handlers/information for job submission, file
  location, and other site-specific matters.

The goal of this tool:

* Develop a research workflow that is general enough and can be run
  across different sites (depending on the need of each calculation,
  where the resources/input files reside for different phase of
  calculations).

* Provide a basic and reliable framework for reproducibility in
  research.
  Researchers write their research procedure as high-level scripts,
  that can be used by others to reproduce their research, or further
  their worn.
  


## Author

Wirawan Purwanto <wirawan0 at gmail.com> .

This library was a byproduct of my research activity at
the College of William and Mary,
funded in part by the National Science Foundation,
Department of Energy, and Office of Naval Research.

WP would like to thank Dr. Henry Krakauer and Dr. Shiwei Zhang at the
College of William and Mary for allowing me to release this library to
the public.

This tool is released with a liberal software license (Apache 2.0) in
order to allow other researchers to pick and reuse the tools that I
have made for their own use.
Contributions back to the toolbox are welcome and strongly encouraged.
