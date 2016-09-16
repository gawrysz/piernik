#!/bin/bash

rm -rf runs/advection_test 2> /dev/null

./setup advection_test -p problem.par.restart_test_v2 

(
    cd runs/advection_test

    echo "run_id = ts1"
    echo "Start:   t = 0.0, nproc = 1"
    mpiexec -n 1 ./piernik > ts1.out
    echo "Restart: t = 1.0, nproc = 1"
    mpiexec -n 1 ./piernik -n '$END_CONTROL tend = 2.0 /' >> ts1.out
    echo "Finish:  t = 2.0"
    echo
    echo "run_is = ts2"
    echo "Start:   t = 0.0, nproc = 5"
    mpiexec -n 5 ./piernik -n '$OUTPUT_CONTROL run_id = "ts2" /' > ts2.out
    echo "Restart: t = 1.0, nproc = 3"
    mpiexec -n 3 ./piernik -n '$END_CONTROL tend = 2.0 /' -n '$OUTPUT_CONTROL run_id = "ts2" /' >> ts2.out
    echo "Finish:  t = 2.0"
)

./bin/gdf_distance runs/advection_test/moving_pulse_ts{1,2}_0002.h5 | tee compare.log

