# silverblaze
### Version 1.0.0
[![Travis build status](https://travis-ci.org/Michael-Stevens-27/silverblaze.svg?branch=master)](https://travis-ci.org/Michael-Stevens-27/silverblaze)
[![AppVeyor build status](https://ci.appveyor.com/api/projects/status/github/Michael-Stevens-27/silverblaze?branch=master&svg=true)](https://ci.appveyor.com/project/Michael-Stevens-27/silverblaze)
[![DOI](https://zenodo.org/badge/127313359.svg)](https://zenodo.org/badge/latestdoi/127313359)

--------------------------------------------------------------------------------------------------------------------------------

Silverblaze offers a finite mixture model for Geographic Profiling. Through the package a user can 

* infer source locations and dispersal patterns - common desirable parameters in geographic profiling
* make inferences via count,  prevalence and point pattern data   
* specify a spatial prior in source locations via shape files

Silverblaze is a toolbox for geographic profiling that utilises different MCMC algorithms (Metropolis-Hastings coupling and Gibbs sampling) to infer the above parameters. The package [RgeoProfile](https://github.com/bobverity/Rgeoprofile) uses a Dirchlet Process Mixture model (DPM) to solve the issue of multiple sources.           
