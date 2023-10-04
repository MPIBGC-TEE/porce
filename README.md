# Photosynthesis, allocation, Organic matter dynamics, and RadioCarbon  Exchange (PORCE)
PORCE is a modelling framework for representing carbon dynamics in terrestrial ecosystems.
It is distributed in the form of an R package, and it uses concepts from functional programming and 
dynamical systems theory, with functions that allow code parallelization for calculation of model uncertainties.
PORCE is built on top of the SOIL-R package, with specific functions for representing photosynthesis, carbon allocation,soil organic matter dynamics, and radiocarbon. 

The functional programming approach allows users to easily exchange different functions for representing specific
processes. An ecosystem model in PORCE can be solved as a steady-state system or as a transient simulation with a given set of initial conditions and forcing variables. The system can also compute solutions as the convolution of an 
impulse response function.  

