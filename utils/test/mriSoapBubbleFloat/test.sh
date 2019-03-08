#!/usr/bin/env bash
source "$(dirname $0)/../../../test.sh"

test_command test_soapbubble src.mgz ctrl.mgz dst.mgz
compare_vol dst.mgz ref.mgz
