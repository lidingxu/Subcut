# Test codes of cutting planes for submodular optimization

We develop cutting planes that could be used for general submodular maximization and submodular-supermodular maximization problems. The cutting planes are integrated into the branch-and-cut algorithm.
The cuts and experiments are described in the paper [Submodular maximization and its generalization through an intersection cut lens](https://arxiv.org/abs/2302.14020). 

The cutting planes are generated in each node of the branch-and-cut algorithm, and we use [SCIP](https://www.scipopt.org/) to implement them. Note that the codes are developed in a Linux environment and the building requires the [Armadillo](https://arma.sourceforge.net/) library. Before running the code, please build them using an appropriate environment.


Three benchmark tests correspond to three root directories. In each benchmark file, there are codes and bash files to run the experiments. 
