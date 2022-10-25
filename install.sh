#!/bin/bash -i

#################################################################
# Install script initialization
#################################################################

# Exit this script if a command fails
set -o errexit

# Constants
RC_FILE="$HOME/.bashrc"
CONDA_ENV_FILE="environment.yml"
REQ_FILE="requirements.txt"

# Check for command-line arguments
while getopts a: flag
do
    case "${flag}" in
        a) aleaf_dir=${OPTARG};;
    esac
done

#################################################################
# Set values for environment variables
#################################################################

# Set ABCE_DIR to the location of this script
abce_dir=$( dirname -- $( readlink -f -- "$0"; ) )
echo "\$ABCE_DIR will be set to $abce_dir"
export ABCE_DIR=$abce_dir

# Set up ALEAF_DIR
# If specified by a command-line argument, don't prompt the user for the
#   ALEAF_DIR location. The empty string "" is a valid command-line value
#   for this directory.
if [[ ! -n ${aleaf_dir+x} ]]; then
    # If not specified in the command line, request a value from the user
    echo "Please enter the absolute path to the top level of the ALEAF source directory."
    echo "If you do not have ALEAF installed, press Enter to leave this variable empty."
    read aleaf_dir
fi
echo "\$ALEAF_DIR will be set to $aleaf_dir"

#################################################################
# Install Julia 1.8, if not already found
#################################################################

# If julia 1.8 is not the current version of julia, download and install
#   julia-1.8.2
updated_julia=0
if [[ -z $( julia --version | grep "1.8" ) ]]; then
    # If julia --version doesn't return something containing the value "1.8",
    #   then julia 1.8 needs to be installed or set to the default local julia
    if [[ ! -d "$HOME/julia-1.8.2" ]]; then
        if [[ ! -f "$HOME/julia-1.8.2-linux-x86_64.tar.gz" ]]; then
            echo "Downloading Julia 1.8.2...";
            wget https://julialang-s3.julialang.org/bin/linux/x64/1.8/julia-1.8.2-linux-x86_64.tar.gz -P "$HOME";
        fi

        echo "Installing Julia 1.8.2..."
        tar zxvf "$HOME/julia-1.8.2-linux-x86_64.tar.gz";
    fi
    updated_julia=1
fi


#################################################################
# Update the .bashrc file
#################################################################

# If the .bashrc file doesn't already exist: create it
if [[ ! -f "${RC_FILE}" ]]; then
    echo "No .bashrc file found; creating a new one at ${RC_FILE}";
    touch "${RC_FILE}";
fi

# Create an associative array to allow looping over ABCE environment variables
declare -A env_vars
env_vars=( ["ABCE_DIR"]="$abce_dir" ["ALEAF_DIR"]="$aleaf_dir" )

# Update or append values for each environment variable in the rc file
echo "Updating environment variables in ${RC_FILE}"
for var_name in "${!env_vars[@]}";
do
    if grep -q "${var_name}" "${RC_FILE}"; then
    # If the form "export $var_name" is found at the start of a line,
    #   replace the rest of the line with the appropriate value
        echo "Found ${var_name}; updating its referenced path";
        sed -i "s|^export $var_name=.*|export $var_name=${env_vars[$var_name]}|" "${RC_FILE}";
    else
    # If the form "export $var_name" does not start any lines in the rc file,
    #   append a line at the end of the rc file to export this variable
        echo "Did not find ${var_name}; adding a new line to .bashrc";
        echo "export ${var_name}=${env_vars[$var_name]}" >> "${RC_FILE}"
    fi
done

# If this script installed Julia 1.8.2, make sure it's added appropriately
#   to the PATH
if [[ $updated_julia ]]; then
    echo "Ensuring Julia 1.8.2 is added to \$PATH in ${RC_FILE}";
    if grep -q "export PATH=.*julia.*" "${RC_FILE}"; then
        # If .bashrc already has a line adding julia to the path, update it to 1.8.2
        echo "Updating existing julia reference in ${RC_FILE}";
        sed -i "s|^export PATH=.*julia-.*|export PATH=\$PATH:$HOME/julia-1.8.2/bin/|" "${RC_FILE}";
    else
        # Otherwise, add a new line adding julia to the path
        echo "Adding julia to \$PATH";
        echo "export PATH=\$PATH:$HOME/julia-1.8.2/bin/" >> "${RC_FILE}";
    fi
fi


#################################################################
# Set up the Python environment
#################################################################

# Determine whether the script is running in a conda environment
# If conda is installed and available for environment management, use it
if [[ ! -z $( conda --version | grep -Eo "conda.*[0-9][0-9]\.[0-9]\.[0-9]" ) && ! -z $( conda info --envs | grep "\*" ) ]]; then
    echo "conda environment detected; using conda to manage python packages"
    CONDA_ENV_FILE="$ABCE_DIR/$CONDA_ENV_FILE"
    # If a conda environment called abce_env doesn't already exist, create it
    if [[ -z $( conda info --envs | grep "abce_env" ) ]]; then
        echo "No conda environment named abce_env found; creating new environment"
        conda env create -f "${CONDA_ENV_FILE}"
    # If this environment already exists, update it
    else
        echo "Found preexisting conda environment named abce_env; updating it"
        conda env update --file "${CONDA_ENV_FILE}" --prune
    fi
# If conda is not available for environment management, use pip to install
#   packages directly
else
    echo "Using pip to manage python packages"
    python3 -m pip install --upgrade pip
    if [[ -f "$ABCE_DIR/$REQ_FILE" ]]; then
        pip install -r "$ABCE_DIR/$REQ_FILE";
    fi
fi

echo "Python environment created successfully."

#################################################################
# Set up the Julia environment
#################################################################


echo "Setting up the local Julia environment..."
julia make_julia_environment.jl

echo "Julia environment created successfully."

#################################################################
# Cleanup
#################################################################

# Prevent current value of $aleaf_dir from contaminating future runs of this
#   script in the same terminal session
unset aleaf_dir
unset abce_dir

echo "All setup completed successfully."
echo "Please close and restart your terminal session so the changes can take effect."
