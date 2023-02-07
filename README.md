# tsp
Models and cut separation callbacks (fractional and integer solutions) for the traveling salesman problem - TSP.

There are two models implemented, i.e. MTZ and the Lazy, which relies on subtour elimination constraints.

CUT and SEC subtour elimination constraints can be chosen through the parameter lazyType in function lazyModel().

lazyModel() also has a parameter named silent which defaults to true. If passed as false, all subtours found and correspondent cuts
added to the formulation are printed.
