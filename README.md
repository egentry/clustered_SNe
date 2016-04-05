# Clustered SNe
Riemann solver with cooling, for evolving a superbubble produced by 1-1000 SNe from a single cluster
-------

Author: Eric Gentry   (gentry.e@gmail.com; egentry@ucsc.edu)   

Licensed under the GPLv3.

-------

##Main Objectives
  - Evolve a spherical symmetric blast wave, incorporating hydrodynamics and radiative cooling
    - Achieve accurate and stable results through Sedov and (thin shell) radiative phases
  - Incorporate multiple supernovae, with realistic delay times and energies
  - Measure the energy and momentum injected into the surrounding medium, to be used as feedback prescriptions in low resolution galaxy/cosmology simulations


##Setup
- Install dependencies listed below
- Adapt `src/makefile` to reflect location of required libraries (`INC_*` and `LIB_*` variables)
- In the `src` directory, run `make all install clean`


##Requires
  - c++ compiler (assumes clang for OS X, gcc for Linux)
    - should (must?) be c++11 capable
  - [slug2](https://bitbucket.org/krumholz/slug2)
    - slug2 must be built as a shared library (change into the `src` directory and call `make lib`; if you enabled FITS at compile time, then keep it enabled at link time. For simplicity, use `make all && make lib` in the `src` directory.)
    - A fork frozen to the version used in my simulations can be found at: https://bitbucket.org/egentry/slug2
    - requires the [GSL](https://www.gnu.org/software/gsl/), [Boost](http://www.boost.org/)
  - libuuid
  - [grackle cooling](https://bitbucket.org/grackle/grackle)
    - requires [HDF5](https://www.hdfgroup.org/HDF5/release/obtain5.html)
  - For visualization:
    - Python (tested for v3.5, maybe backwards compatible)
    - Jupyter notebook (tested for v4)
      - ipywidgets (tested for v4)
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


##Acknowledgements
This project built upon [RT1D](https://github.com/duffell/RT1D), an open-source riemann solver from [Paul Duffell](http://duffell.org/).

This project was also includes the `sedov3.f` code of http://cococubed.asu.edu/research_pages/sedov.shtml in order to generate the analytic sedov solution for the purpose of visualization.  While it's not difficult to generate a basic sedov solution, it's difficult to do well. Using `sedov3.f` allows us to avoid worrying about the details of quad precision math ourselves.
