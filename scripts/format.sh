#!/usr/bin/env bash
# Format src code using the pep8 style guide.
SCRIPTS_DIR="$(dirname $(realpath $0))"
REPO_DIR="$(dirname $SCRIPTS_DIR)"
PACKAGE_NAME="simdata"
SRC_DIR="$REPO_DIR/src/$PACKAGE_NAME"
TEST_DIR="$REPO_DIR/tests"
if [[ -d "$SRC_DIR" ]]; then
	yapf --recursive --in-place --style pep8 "$SRC_DIR"
fi
if [[ -d "$TEST_DIR" ]]; then
	yapf --recursive --in-place --style pep8 "$TEST_DIR"
fi
