# Cmake version needed to support -B -S
# cmake .. -DSCIP_DIR=/home/ldx/Projects/scip-8TEST -DARMA_DIR=/home/ldx/Projects/armadillo
#!/bin/bash
SCIP_DIR=/home/ldx/Projects/scip-8TEST
ARMA_DIR=/home/ldx/Projects/armadillo

# install D optimal
rm -rf  Doptimal/solver/build/*


echo 
# install max cut
rm -rf  Maxcut/solver/build/*


echo

# install multilinear
rm -rf  Multilinear/solver/build/*
