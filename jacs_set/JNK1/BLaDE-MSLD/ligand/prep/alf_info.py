import shutil
import subprocess
import numpy as np
import os
import re

def get_executable_or_alias(name):
    # Check if it's in $PATH
    path = shutil.which(name)
    if path:
        return path

    # If it's not in $PATH, check if it's an alias
    result = subprocess.run(['bash', '-i', '-c', f'alias "{name}"'], capture_output=True, text=True)
    if result.stdout:
        alias = result.stdout
        #return alias
        s=alias
        match = re.search(r"alias blade=\'(?P<path>.+?)\'", s)
        if match:
           path = match.group('path')
           print(path)
           return path   

    return None

# usage
path_or_alias = get_executable_or_alias('blade')

blade_env = os.getenv('blade')
if blade_env != None: 
   path_or_alias = blade_env


alf_info={}
alf_info['name']='system'
alf_info['ncentral']=0
alf_info['nreps']=1
alf_info['nnodes']=1
alf_info['enginepath']=path_or_alias
alf_info['temp']=300.0
alf_info['nsubs']=[3, 4, 4, 4]
alf_info['nblocks']=np.sum(alf_info['nsubs'])
