#!/bin/bash

# Use testSysimage.so if available
if [[ -f "$ABCE_DIR/test/testSysimage.so" ]]; then
    julia --project="$ABCE_DIR/env" -J "$ABCE_DIR/test/testSysimage.so" "$ABCE_DIR/test/test_julia.jl"
else
    julia --project="$ABCE_DIR/env" "$ABCE_DIR/test/test_julia.jl"
fi
