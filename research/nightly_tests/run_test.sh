#!/bin/bash

# This code runs a moris package test and checks its output

# Input arguments
NIGHTLY_TEST_DIR=$1
PACKAGE=$2

# Set variables
WORKING_DIRECTORY=$NIGHTLY_TEST_DIR/$PACKAGE/build

# Get the system enivironment variables
# source /home/ryan/.bashrc

# Move to working directory
cd $WORKING_DIRECTORY

# Run tests
ctest -VV >&1
ctest --rerun-failed -VV

# Prevent errors throwing a stop to parent script
exit 0
