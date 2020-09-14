#!/bin/bash

prime=$1
algorithm=$2

( python3 parameters.py -p $prime -f svelu -a $algorithm ) | tee tmp/$prime

( tail -n +3 tmp/$prime ) | tee ijk/$prime

rm tmp/$prime
