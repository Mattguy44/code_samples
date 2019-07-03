#!/bin/bash

# The following script builds a cmake package and tests if it built without errors.
# It is assumed that this script has been called from build_test_runner.py.

# Input arguments
NIGHTLY_TEST_DIR=$1
PACKAGE=$2
CMAKE_BUILD_FLAGS="$3 $4"
echo "Flags: $CMAKE_BUILD_FLAGS"
# Set variables
WORKING_DIR=$NIGHTLY_TEST_DIR/$PACKAGE

# Get the system enivironment variables
# source /home/ryan/.bashrc

# Move to package working directory
cd $WORKING_DIR

# Remove old package build
rm -rf build/
mkdir build

# Move to build directory
cd build

# Build package
cmake $WORKING_DIR/$PACKAGE $CMAKE_BUILD_FLAGS > $WORKING_DIR/logs/cmake_output.txt
cmake $WORKING_DIR/$PACKAGE # resolves cashed variables
make -j6 >& $WORKING_DIR/logs/make_output.txt 2> $WORKING_DIR/logs/err_output.txt

# Find any errors in the build
grep "error:" "$WORKING_DIR"/logs/err_output.txt >> "$WORKING_DIR"/logs/shortFails.out

# If grep found "error:", then $? will return 0.
# The command $? returns 0 if the previous command was successful and 1 otherwise.
BUILD_FAILED=$?

# Return an error message on failure
echo "@@~~@@"
if [[ $BUILD_FAILED -eq 0 ]]
then
    echo "The build encountered errors with the following commands:"
    echo "cmake $WORKING_DIR/$PACKAGE $CMAKE_BUILD_FLAGS"
    echo "make -j6"
else
    echo "The build succeeded with the following commands:"
    echo "cmake $WORKING_DIR/$PACKAGE $CMAKE_BUILD_FLAGS"
    echo "make -j6"
fi
