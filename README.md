# Test codes for cutting planes for signomial programming

We develop two kinds of cutting planes for signomial  programming: outer approximation cuts and intersection cuts.
The cuts and experiments are described in the paper [Submodular maximization and its generalization through an intersection cut lens](https://arxiv.org/abs/2302.14020). 

The cutting planes are generated in each node of the branch and bound algorithm, and we use [SCIP](https://www.scipopt.org/) to implement them. Note that the codes are developed in a Linux environment and the building requires the [Armadillo](https://arma.sourceforge.net/) library.

There are three benchmarks. In each benchmark file, there are codes and bash files to run the experiments. Before running the code, please build them using an appropriate environment.
