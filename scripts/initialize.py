import os, sys, argparse, glob, subprocess

if __name__ == "__main__":
    
    args = sys.argv
    systemPath = os.path.dirname(os.path.abspath(__file__))
    
    parser = argparse.ArgumentParser()
    parser.add_argument('call')
    parser.add_argument('--clean', action='store_true', help='Deletes processed data and outputs prior to run')
    
    cmdargs = parser.parse_args(args)
    if cmdargs.clean:
        for f in glob.glob(os.path.join(systemPath,'..','figures','*.png')):
            os.remove(f)
        for f in glob.glob(os.path.join(systemPath,'..','figures','*.csv')):
            os.remove(f)
        for f in glob.glob(os.path.join(systemPath,'..','output','*.csv')):
            os.remove(f)
        for f in glob.glob(os.path.join(systemPath,'..','data','*_processed','*.csv')):
            os.remove(f)
    
    subprocess.call([sys.executable, os.path.join(systemPath, 'fsri_collect_thermophysical_properties.py')])
    subprocess.call([sys.executable, os.path.join(systemPath, 'fsri_collect_cone_data.py')])
    subprocess.call([sys.executable, os.path.join(systemPath, 'process_fsri_database.py')])
    subprocess.call([sys.executable, os.path.join(systemPath, 'process_faa_data.py')])
    subprocess.call([sys.executable, os.path.join(systemPath, 'process_fpl_data.py')])
    