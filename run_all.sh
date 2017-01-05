#!/bin/bash

for D in `ls jobs`
do
    cd "jobs/$D"
    ../../run_job.R
    cd ../..
done
