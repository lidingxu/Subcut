# Cmake version needed to support -B -S
# cmake .. -DSCIP_DIR=/home/ldx/Projects/scip-8TEST -DARMA_DIR=/home/ldx/Projects/armadillo
#!/bin/bash
SCIP_DIR=/home/ldx/Projects/scip-8TEST
ARMA_DIR=/home/ldx/Projects/armadillo

# install D optimal
rm -rf  Doptimal/solver/build/*
cmake -B Doptimal/solver/build -S Doptimal/solver -DSCIP_DIR=$SCIP_DIR -DARMA_DIR=$ARMA_DIR
make -C Doptimal/solver/build

echo 
# install max cut
rm -rf  Maxcut/solver/build/*
cmake -B Maxcut/solver/build -S Maxcut/solver -DSCIP_DIR=$SCIP_DIR
make -C Maxcut/solver/build

echo

# install multilinear
rm -rf  Multilinear/solver/build/*
cmake -B Multilinear/solver/build -S Multilinear/solver -DSCIP_DIR=$SCIP_DIR
make -C Multilinear/solver/build
