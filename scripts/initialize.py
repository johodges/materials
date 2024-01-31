import os, sys, argparse, glob, subprocess

if __name__ == "__main__":
    
    args = sys.argv
    systemPath = os.path.dirname(os.path.abspath(__file__))
    
    parser = argparse.ArgumentParser()
    parser.add_argument('call')
    parser.add_argument('--clean', action='store_true', help='Deletes processed data and outputs prior to run')
    parser.add_argument('--donotinitialize', action='store_true', help='Deletes processed data and outputs prior to run')
    
    cmdargs = parser.parse_args(args)
    if cmdargs.clean:
        for f in glob.glob(os.path.join(systemPath,'..','input_files','*','*')):
            os.remove(f)
        dirs = sorted(glob.glob(os.path.join(systemPath,'..','input_files','*')))
        for f in dirs:
            os.rmdir(f)
        for f in glob.glob(os.path.join(systemPath,'..','figures','*.png')):
            os.remove(f)
        for f in glob.glob(os.path.join(systemPath,'..','figures','*.xlsx')):
            os.remove(f)
        for f in glob.glob(os.path.join(systemPath,'..','output','*.csv')):
            os.remove(f)
        for f in glob.glob(os.path.join(systemPath,'..','data','*_processed','*.csv')):
            os.remove(f)
        for f in glob.glob(os.path.join(systemPath,'..','data','*_processed','*','*.csv')):
            os.remove(f)
    
    if cmdargs.donotinitialize is not True:
        my_env = os.environ.copy()
        #subprocess.call([sys.executable, os.path.join(systemPath, 'fsri_collect_thermophysical_properties.py')])
        #subprocess.call([sys.executable, os.path.join(systemPath, 'fsri_collect_cone_data.py')])
        #subprocess.call([sys.executable, os.path.join(systemPath, 'process_fsri_database.py')])
        #subprocess.call([sys.executable, os.path.join(systemPath, 'process_faa_data.py')])
        #subprocess.call([sys.executable, os.path.join(systemPath, 'process_fpl_data.py')])
        print("Step 1/5: Collecting fsri thermophysical material properties")
        process = subprocess.Popen([sys.executable, os.path.join(systemPath, 'fsri_collect_thermophysical_properties.py')], env=my_env, shell=False, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        out, err = process.communicate()
        errcode = process.returncode  
        print("Step 2/5: Collecting fsri cone data")
        process = subprocess.Popen([sys.executable, os.path.join(systemPath, 'fsri_collect_cone_data.py')], env=my_env, shell=False, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        out, err = process.communicate()
        errcode = process.returncode  
        print("Step 3/5: Processing fsri data")
        process = subprocess.Popen([sys.executable, os.path.join(systemPath, 'process_fsri_database.py')], env=my_env, shell=False, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        out, err = process.communicate()
        errcode = process.returncode  
        print("Step 4/5: Processing faa data")
        process = subprocess.Popen([sys.executable, os.path.join(systemPath, 'process_faa_data.py')], env=my_env, shell=False, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        out, err = process.communicate()
        errcode = process.returncode  
        print("Step 5/5: Processing fpl data")
        process = subprocess.Popen([sys.executable, os.path.join(systemPath, 'process_fpl_data.py')], env=my_env, shell=False, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        out, err = process.communicate()
        errcode = process.returncode  
        print("Initialization complete")
    