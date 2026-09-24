#!/bin/sh
# Validation suite for the split census (docs/CN2_Exhaustive_Note.md §9.3).
#   tests/run_validation.sh          quick suite (a few minutes)
#   tests/run_validation.sh --full   adds the heavier checks (about 45 minutes more on 30 threads; the brute force needs 40 GB RAM)
cd "$(dirname "$0")/.." && exec python3 tests/validate.py "$@"
