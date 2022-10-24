#!/bin/bash

# Exit this script if a command fails
set -o errexit

# Check for command-line arguments
while getopts a: flag
do
    case "${flag}" in
        a) aleaf_dir=${OPTARG};;
    esac
done

echo $aleaf_dir

# Set up environment variables
# Set ABCE_DIR to the location of this script
ABCE_DIR=$( dirname -- $( readlink -f -- "$0"; ) )
echo "Setting environment variable ABCE_DIR=$ABCE_DIR"
export ABCE_DIR=$ABCE_DIR

# Set up ALEAF_DIR
# If specified by a command-line argument, don't prompt the user for the
#   ALEAF_DIR location
if [[ ! -n ${aleaf_dir+x} ]]; then
    # If not specified in the command line, request a value from the user
    echo "Please enter the absolute path to the top level of the ALEAF source directory."
    echo "If you do not have ALEAF installed, press Enter to leave this variable empty."
    read aleaf_dir
fi

echo "Setting environment variable ALEAF_DIR=$aleaf_dir"
export ALEAF_DIR=$aleaf_dir


# Set up the Julia environment
echo "Setting up the local Julia environment..."
julia make_julia_environment.jl --clean

echo "Julia environment set up."

# Cleanup
unset aleaf_dir

echo "All setup completed successfully."
