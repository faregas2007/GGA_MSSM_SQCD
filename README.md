# GGA_MSSM_SQCD
Individual run - by default will use python3, if one use cmd mode 

cd input_dir

python3 {core_dir}/run.py -id {diagrams_id} -lam {lambda_shift} -in {input_file} -amp {amp_dir}

The json file can be configurated in the way that you want. In order to use Richardson extrapolation, richardson flag needs to be set to true. Then, one can change the k, order of the Richardson extrapolation to the expected order. By default, The calculation will need a cluster to take the advatanges of parallel computing. Using 

python3 main.py

All of the amplitudes here are numerical stable below and above quark threshold upto lam=10^-4,-5. The amplitudes of diagram 2, 12, 13 and 14 will become unstable above squark threshold. 
