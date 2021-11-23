#!/bin/sh
#
# A pre-push test script for ABCE
# Tests to see whether an ABCE run crashes with current code
#   configuration and settings.

# Describe the test to the user
echo "TEST: ABCE simple crash test"
echo "This test will run ABCE using whatever settings currently exist in settings.yml."
echo "If ABCE does not crash (i.e. it runs to completion and returns an exit code of 0, the test passes."
echo "If something happens and ABCE returns any non-zero return code, the test fails.\n"
echo "Warning: This test does NOT guarantee code is working correctly! It only detects fatal crashes."

# Use the CL argument to set the path to run.py
# If no CL was given, throw an error and exit with return code 1
if [ -z "${1+x}" ]
then
    run_script_path="$1"
else
    echo "No path specified for the crash test; aborting."
    exit 1
fi

# Run the test
run_cmd="python3 $run_script_path -f"
${run_cmd}

