#!/usr/bin/env bash

set -euox pipefail

bioino --help

bioino gff2table test/test.gff3 > /dev/null
