# silverblaze
### Version 0.0.1 (pre-release)
[![Travis build status](https://travis-ci.org/Michael-Stevens-27/silverblaze.svg?branch=master)](https://travis-ci.org/Michael-Stevens-27/silverblaze)
[![AppVeyor build status](https://ci.appveyor.com/api/projects/status/github/Michael-Stevens-27/silverblaze?branch=master&svg=true)](https://ci.appveyor.com/project/Michael-Stevens-27/silverblaze)

--------------------------------------------------------------------------------------------------------------------------------

The R package silverblaze builds on the foundations of the geographic profiling methododlogy by implementing information drawn from absence data. Using the locations of sentinel sites, each associated with a capture density, the model infers:

1. The most likely areas that contain a pre-specified number of sources.
2. The dispersal parameter sigma.
3. The underlying population density.

This presence-absence model utilises MCMC algorithms (Metropolis-Hastings and Gibbs sampling) to infer these parameters. As stated silverblaze builds on an existing methododlogy. The package [RgeoProfile](https://github.com/bobverity/Rgeoprofile) uses a Dirchlet Process Mixture model (DPM) to solve the issue of multiple sources.           
