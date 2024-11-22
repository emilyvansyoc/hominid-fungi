#!/bin/bash

##### run Andrew's script to generate random trees

# activate conda environment
source ~/miniconda3/etc/profile.d/conda.sh
conda activate py2

## ---- set variables ----

# set working directory
DIR=/mypath/

# number of random trees to generate
NRAND=100000

# tree labels (comma separated list with no spaces)
LABS=Gorilla_beringei,Gorilla_gorilla,Homo_sapien,P_troglodytes_schweinfurthii,Mangabey

# name for random tree (outside variable to be easy to run)
RAND=$1

## ---- run random trees ----

# make temporary directory for random trees
TMP=$DIR/TEMP
mkdir -p $TMP

# run random trees
python /mypath/tree_random.py -n $NRAND -i $LABS -o $TMP

# concatenate trees into one file
for i in $TMP/*
do echo $i
    cat $i >> $RAND
done

# remove temp directory
rm -rf $TMP

# print finished message
echo "done making random trees"