mpirun.mpich -np 6 python -W ignore 0_preprocessing.py > 0_preprocessing.txt
echo $"\n\nRamping up...\n"
mpirun.mpich -np 6 python -W ignore 1_ramp.py > 1_ramp.txt
echo $"\n\nRunning...\n"
mpirun.mpich -np 6 python -W ignore 2_run.py > 2_run.txt
