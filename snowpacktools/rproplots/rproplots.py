
import os
import subprocess
import multiprocessing as mp
import pkg_resources

from datetime import datetime

def plot_single_profile_R(pro_file, outdir, dt_nowcast_end, daily_time="23:00", dt_forecast_end=None, smet_file=None, i=0, **kwargs):

    os.makedirs(outdir, exist_ok=True)
    datetime_format   = '%Y-%m-%dT%H:%M'
    dt_nowcast_end = datetime.strftime(datetime.strptime(dt_nowcast_end, datetime_format), datetime_format)
    command = ["Rscript", pkg_resources.resource_filename('snowpacktools', 'rproplots/plot_profiles.R'), 
               "--pro", pro_file,
               "--outdir", outdir,
               "--dt_nowcast_end", str(dt_nowcast_end),
               "--daily_time", daily_time
               ]
    
    if dt_forecast_end:
        dt_forecast_end = datetime.strptime(dt_forecast_end, datetime_format)
        command.extend(["--dt_forecast_end", str(dt_forecast_end)])
    if smet_file:
        command.extend(["--smet_file", smet_file])

    try:
        result = subprocess.run(
            command, capture_output=True, text=True, check=True
        )
        return result.returncode
    except subprocess.CalledProcessError as e:
        print(
            f"[E]  R plot script failed for {os.path.basename(pro_file)} with return code {e.returncode}.\n"
            f"     Command: {' '.join(e.cmd)}\n"
            f"     Output: {e.stdout}\n"
            f"     Error: {e.stderr}"
        )
        return e.returncode


def plot_profiles_R(pro_files, outdir, dt_nowcast_end, ncpus=mp.cpu_count()-2, **kwargs):

    args_list = [(pro, outdir, dt_nowcast_end, *kwargs) for pro in pro_files]
    with mp.Pool(processes=ncpus) as pool:
        results = pool.starmap(plot_single_profile_R, args_list)
    
    failed = [pro_files[i] for i, return_code in enumerate(results) if return_code != 0]
    if failed:
        print(f"[E]  R plot script failed for the following profiles: {failed}")
    else:
        print(f"Plotted all profiles successfully in {outdir}.")
    
