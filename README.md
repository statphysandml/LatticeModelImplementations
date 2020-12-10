LatticeModelImplementations
=================

LatticeModelImplementations is C++ library which is build upon the MCMCSimulationLib. It implements different well-known Monte Carlo algorithms and interesting models in statistical physics. In the long run, the library should be a collection for all common used models. Switching between different models and different algorithms is due to the modular structure very easy.

The evaluation of the simulation data takes place in python. The python modules allow a convenient computation of observables and further interesting quantities. In addition, the python modules can be used to transform the data into a pytorch dataset and to train a machine learning algorithm on the data.