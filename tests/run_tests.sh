#!/bin/bash
# cofactor's own regression tests: extractor and CLI (Python), check/report library (C++, plain g++ or clang++).
# No Enzyme needed. Exit 0 only if everything passes.
set -u
cd "$(dirname "$0")/.."
CXX=${CXX:-g++}
TMP=$(mktemp -d)
rc=0
echo "== python tests"
python3 -m unittest discover -s tests -p 'test_*.py' -v || rc=1
echo "== check/report library"
$CXX -std=c++14 -Wall -Ilib tests/test_checks.cc -o $TMP/t_checks && $TMP/t_checks || rc=1
echo "== check/report library, forward-only build (incomplete runs exit 2)"
$CXX -std=c++14 -Wall -Ilib -DFWD_ONLY tests/test_checks.cc -o $TMP/t_checks_fwd && $TMP/t_checks_fwd || rc=1
rm -rf $TMP
echo "== $([ $rc = 0 ] && echo ALL TESTS PASSED || echo TESTS FAILED)"
exit $rc
