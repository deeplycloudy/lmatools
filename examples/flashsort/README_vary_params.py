
#from within a working, case-specific subdirectory of the directory containing autorun,
export PYTHONPATH=$PYTHONPATH:`cd ..; pwd`

#Edit vary_params.py
python vary_params.py
#Edit path in examine_big_run.py
python examine_big_run.py