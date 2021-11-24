#!/bin/sh
#
# A pre-push test script for ABCE
# Tests to see whether an ABCE run crashes with current code
#   configuration and settings.

# Describe the test to the user
printf "TEST: ABCE simple crash test\n\n"
printf "This test will run ABCE using whatever settings currently exist in settings.yml.\n"
printf "If ABCE does not crash (i.e. it runs to completion and returns an exit code of 0, the test passes.\n"
printf "If something happens and ABCE returns any non-zero return code, the test fails.\n\n"
printf "Warning: This test does NOT guarantee code is working correctly! It only detects fatal crashes.\n\n"

# Set pre-push test settings file name
test_settings_file="$ABCE_DIR/pre_push_test_settings.yml"

# Make a test-specific copy of settings.yml
cp "$ABCE_DIR/settings.yml" $test_settings_file

# Change the `num_steps` value to 2
sed -i 's|num_steps:.*|num_steps: 2|' $test_settings_file

# Run the test
run_cmd="python3 $ABCE_DIR/run.py -f --settings_file=$test_settings_file"
${run_cmd}

# Save ABCE's return code while we do something else
ret_code=$?

# Delete the test copy of settings.yml
rm "$test_settings_file"

exit $ret_code
