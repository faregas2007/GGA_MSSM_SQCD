# each qcd part will render the sqcd part. Since the qcd part only depends on quarks and higgs masses.
# submit only sqcd parts
# for multiple inputfiles, only qcd parts are needed to submit to cluste 

# for now, do it manually. Generate input--> csqcd--> xs 

def submit(
    jobname: str: params.sushi_input,
    men : int=8,
    inputfile: str=Path(working_dir,params.sushi_input)
    outputfile: str=Path(working_dir, params.sushi_input)
    src_dir: str=src_dir,
    sqcd: bool=False
):
    if(sqcd==False):
        wrap = f"python3 {src_dir}/run.py"
    else:
        wrap = f"python3 {src_dir}/sqcd_run.py -id {index_in} -lam {lam} -in {inputfile} -amp {amp_dir}"

    logger.info(f'currently on f<cluster::submit> {wrap}')
    sub = [
        'sbatch',
        '-p albatros',
        f'--job-name={jobname}',
        f'--output="{output_file}"',
        f'--mem="{mem}"',
        f'--wrap="{wrap.strip()}'
    ]

    process = sps.Popen("".join(sub), shell=True, stdout=sps.PIPE)
    stdout = process.communicate()[0].decode('utf-8')

def cluster_run(
    inputfile: str=params.inputfile,
    outputfile: str=params.outputfile,
    src_dir: str=src_dir,
    filename: str=params.filename
    ):
    # the only difference would be inputfile_name
    if(sqcd == True):
        jobname = filename.strip('.csv')
        submit(
            jobname =jobname,
            mem=8,
            inputfile=inputfile,
            outputfile=outputfile,
            src_dir=src_dir
        )
        logger.info(f'Job {jobname} is submitted.')
    else:
        if(Micheal == True):
            jobname = 'xs_'+filename.strip('.in')
        else:
            jobname = 'xs_'+filename
        submit(
            jobname=jobname,
            mem=8,
            inputfile=inputfile,
            outputfile=outputfile,
            src_dir=src_dir
        )
        logger.info(f'Job {jobname} is submitted.')
    return

if __name__ == "__main__":
    jobs = ['test_xs1_nonbatch', 'test_xs2_nonbatch']
    for job in jobs:
        submit(
            jobname=job,
            inputfile = params.sushi_input,
            outputfile = Path(working_dir, 'xs_test'+params.sushi_input),
            sqcd = False
        )
    
    """
    if os.path.isfile(filenamer)
        csvwriter(
            working=str(working_dir)
            name = params.name
        )
        break
    else:
        time.sleep(60)
    
    logger.info(f'xs evaluation on {params.name} is done.')
    logger.info(f'The result is stored at {working_dir}/{filename}')
    """
