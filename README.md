# Test codes for cutting planes for signomial programming

The cuts and experiments are described in the paper [Submodular maximization and its generalization through an intersection cut lens](https://arxiv.org/abs/2302.14020). We test two kinds of cuts: outer approximation cuts and intersection cuts.


The cutting planes are generated in each node of the branch and bound algorithm, and we use [SCIP](https://www.scipopt.org/) to implement it.

There are three benchmarks. In each benchmark file, there is a code and bash file to run the experiments. Before running the code, please build them using an appropriate environment.
