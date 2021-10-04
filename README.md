# GGA_MSSM_SQCD
Individual run - by default will use python3, if one use cmd mode 

cd input_dir

python3 {core_dir}/run.py -id {diagrams_id} -lam {lambda_shift} -in {input_file} -amp {amp_dir}

Change the json file to the configuration that you want. You have turn on richardson = true, and change the k, order of the Richardson extrapolation, to whatever order that you want. By default, it will need a cluster to take advatanges parallel computing instead of individual run. Using 

python3 main.py

All of the amplitudes are numerical stable below and above quark threshold upto lam=10^-4,-5. These amplitudes will become unstable above squark threshold, for whatever values of lam. 