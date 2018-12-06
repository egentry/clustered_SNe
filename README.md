# Clustered SNe
Riemann solver with cooling, for evolving a superbubble produced by 1-1000 SNe from a single cluster
-------

Author: Eric Gentry   (gentry.e@gmail.com; egentry@ucsc.edu)   

Licensed under the GPLv3.

-------

## Main Objectives
  - Evolve a spherical symmetric blast wave, incorporating hydrodynamics and radiative cooling
    - Achieve accurate and stable results through Sedov and (thin shell) radiative phases
  - Incorporate multiple supernovae, with realistic delay times and energies
    - Include pre-SNe winds as a constant wind through a star's life
  - Measure the energy and momentum injected into the surrounding medium, to be used as feedback prescriptions in low resolution galaxy/cosmology simulations


## Getting started
- Install dependencies listed below
- Adapt `src/makefile` to reflect location of required libraries (`INC_*` and `LIB_*` variables)
- In the `src` directory, run `make all install clean`

That should get you a working executable. Once you have data to process, this repo's wiki has [instructions](https://github.com/egentry/clustered_SNe/wiki/Getting-Starting-with-the-Analysis-Package) for using the python analysis code.


## Requires
  - c++ compiler (assumes clang for OS X, gcc for Linux)
    - should (must?) be c++11 capable
    - you must use the same compiler as you used for installing Boost. If you're unsure, then use the default I set.
  - [slug2](https://bitbucket.org/krumholz/slug2)
    - slug2 must be built as a shared library (change into the `src` directory and call `make lib`; if you enabled FITS at compile time, then keep it enabled at link time. For simplicity, use `make all && make lib` in the `src` directory.)
    - A fork frozen to the version used in my simulations can be found at: https://bitbucket.org/egentry/slug2
    - requires the [GSL](https://www.gnu.org/software/gsl/), [Boost](http://www.boost.org/)
  - libuuid
  - [grackle cooling](https://bitbucket.org/grackle/grackle) (v3)
    - requires [HDF5](https://www.hdfgroup.org/HDF5/release/obtain5.html)
  - For visualization:
    - Python (tested for v3.5, mostly backwards compatible)
    - Jupyter notebook (tested for v4)
      - ipywidgets (tested for v4, v5 -- mildly broken on v6)
    - Matplotlib
    - Bokeh
    - Numpy
    - Scipy
    - Astropy
    - Pandas
    - Seaborn
    - Numba
    - SQLAlchemy
    - corner
      - Only needed if you call `BayesianFit.create_corner_plots`


## Acknowledgements
This project built upon [RT1D](https://github.com/duffell/RT1D), an open-source riemann solver from [Paul Duffell](http://duffell.org/).

This project was also includes the `sedov3.f` code of http://cococubed.asu.edu/research_pages/sedov.shtml in order to generate the analytic sedov solution for the purpose of visualization.  While it's not difficult to generate a basic sedov solution, it's difficult to do well. Using `sedov3.f` allows us to avoid worrying about the details of quad precision math ourselves.

Ben Wibking was a big help in proposing efficiency improvements and identifying bugs.
