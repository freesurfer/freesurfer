#! /bin/sh

gprof -q scuba gmon.out > profile-call.txt
gprof scuba gmon.out > profile-flat.txt
