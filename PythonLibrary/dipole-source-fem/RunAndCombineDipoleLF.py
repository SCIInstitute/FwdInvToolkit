import os
import sys, getopt
import glob

import numpy as np
import scipy.io
import subprocess

thresh = 200


# SCIRun needs to be built with the matrix fix to work
default_scirun = "/Applications/SCIRun.app/Contents/MacOS/SCIRun"

default_tag = "pointsource"
default_volmesh="tetvolmesh_3.59E+03.mat"
default_surfmesh="trisurfmesh_3.59E+03.mat"
default_sourcepoints="volume_source_points_3.59E+03.mat"
default_net_path = "../../Networks/dipole-source-fem"
default_output_path="../../Data/DipoleSphere"

#surfmesh="trisurfmesh"+volmesh[10:]
#pointvol="volume_source_points"+volmesh[10:]

#tmp=scipy.io.loadmat(os.path.join(output_path, pointvol))
#pv_size = tmp["SCIRUNFIELD"]["node"].shape

#tmp=scipy.io.loadmat(os.path.join(output_path, surfmesh))
#tv_size = tmp["SCIRUNFIELD"]["node"].shape

#pv_size=(6642, 1)
#ts_size=(4080, 1)

#pv_size=(13264, 1)
#ts_size=(1296, 1)

pv_size=(1028, 1)
ts_size=(642, 1)

param_dict = {}
param_dict["pointsource"] = {
  "sr_net" : "LLrecipricol_single.srn5",
  "tag"    : "_recip_pointsource_lead",
  "fnum"   : ts_size[0],
  "ssize"  : (pv_size[0], 3)
}
param_dict["forward"] = {
  "sr_net" : "LLdirect_single.srn5",
  "tag"    : "_forward_3vect",
  "fnum"   : pv_size[0],
  "ssize"  : (ts_size[0], 3)
}



def check_sol_size(sol, expectsize):
  return sol.shape == expectsize
  
def check_sol(sol, expectsize, prev_sol = []):
  message ="no problem detected"
  if not prev_sol.any():
    prev_sol = np.zeros(sol.shape)
  
  
  if abs(np.sum(sol)) < 0.000001:
    message = "low amplitude (zeros) detected"
    return False, message
  elif np.max(np.abs(sol))> 1e300:
    message = "possible inf detected"
    return False, message
    
  rms = np.sqrt(np.sum(sol*sol))/np.prod(expectsize)
  ratio = abs(np.std(sol)/rms)
#  ratio = abs(np.std(sol)/np.mean(sol))
#  print("sum: ", np.sum(sol))
#  print("max: ", np.max(np.abs(sol)))
#  print("std: ", np.std(sol))
#  print("mean: ", np.mean(sol))
#  print("RMS: ", rms)
#  print("ratio: ", ratio)
    
  if ratio>thresh:
    message = "possible problem, high std/mean.  ratio: "+str(ratio)
    return False, message
    
  if (sol == prev_sol).all():
    message = "solution same as prev"
#    print(prev_sol)
    return False, message
  
  return True, message
  
  
  
def check_files(sol_path, file_root, num_files):
  files = sorted(glob.glob(os.path.join(sol_path, "*"+file_root+"*")))
  
#  if len(files) == 0:
#    print(" input directory does not have any valid files ('*"+file_root+"*')")
#    return files

  expect_ind = list(range(num_files))
  missing_ind = []
  found_ind = []
  extra_ind = []
  extra_files = []
  cleared_files = []
  
  if len(files) > num_files:
    print("Something went wrong.  Expected " + str(num_files)+ ", but got more " + str(len(files)) + "('*"+file_root+"*')")
  
  if len(files) < num_files:
    print("Missing some files, there should be  " + str(num_files)+ ", but there are " + str(len(files)) + "('*"+file_root+"*')")
      
  for f_full in files:
    
    f_ind = find_index_from_filename(f_full, file_root)
    
    if not type(f_ind)==int:
      continue
    
    found_ind.append(f_ind)
    if not f_ind in expect_ind:
      extra_ind.append(f_ind)
      extra_files.append(f_full)
    else:
      cleared_files.append(f_full)
    
  missing_ind = [ ind for ind in expect_ind if not ind in found_ind ]
  
  if len(missing_ind)>0:
    print("missing indices:")
    print(missing_ind)
  if len(extra_ind)>0:
    print("extra files with same file root provided")
    print(extra_files)
  missing_files = []
  for m_i in missing_ind:
    num_digits = len(str(num_files))
    m_file= file_root+"_"+f"{m_i:05}"+".mat"
#    print(m_file)
    missing_files.append(os.path.join(sol_path, m_file))
    
  return cleared_files, missing_files
  
def load_solutions(files, tag,  expectsize):
  if tag=="forward":
    LF_mat=np.zeros((expectsize[0], expectsize[1]*len(files)))
  else:
    LF_mat=np.zeros((len(files), expectsize[0]*expectsize[1]))
    
  print(LF_mat.shape)
  
  prev_sol = np.zeros(expectsize)
  rerun_files = []
  
  for k, f in enumerate(files):
    tmp = scipy.io.loadmat(f)
    vars=[k for k in tmp.keys() if not k[0]=='_']
    if not vars:
        print("No data found in file",os.path.join(f))
        return False
    elif len(vars)>1:
        print("multiple variables found in file, using: ",vars[0])
        
    sol = tmp[vars[0]]
#    print(sol.shape)
    
    if check_sol_size(sol, expectsize):
      solution_check = check_sol(sol, expectsize, prev_sol)
      if not solution_check[0]:
        print("problem detected with "+f+".  "+solution_check[1])
        rerun_files.append(f)
      
#      print(sol)
      
      if tag=="forward":
        LF_mat[:,k*expectsize[1]:(k+1)*expectsize[1]] = sol
      else:
        LF_mat[k,:] = sol.flatten()
    else:
      print("Something wrong with the size of"+f+".  continuing, but don't use LF")
      rerun_files.append(f)
      prev_sol = np.zeros(expectsize)
      
    prev_sol = sol
  print(LF_mat.shape)
  
  return LF_mat, rerun_files

def find_index_from_filename(filename, file_root):

  f = os.path.basename(filename)
  ind_search = f.split(file_root)
#    print(ind_search)
  
  f_ind = []
  
  for ids in ind_search:
    newstr = ''.join((ch if ch in '0123456789' else ' ') for ch in ids)
    for str_split in newstr.split():
      f_ind.append(int(str_split))
    
#    print(f_ind)
  if len(f_ind)>1:
    print("There may be a problem with the file name: "+f+". Found multiple indices, taking the first")
  elif len(f_ind) == 0:
    print("There may be a problem with the file name: "+f+". No indices found.  Skipping")
    return f_ind
  
  return f_ind[0]

def rerun_solution(filename,file_root,scirun_net, scirun_call, interactive, tvmesh_file, tsmesh_file, vpmesh_file):
  
  f_ind = find_index_from_filename(filename, file_root)
  
  f_path =os.path.dirname(filename)
  ind_file = os.path.join(f_path, file_root+"index_file.txt")
  
  with open(ind_file,'w') as fid:
    fid.write(str(f_ind))
  fid.close()
  
  os.environ["INDEXFILE"]=ind_file
  os.environ["VOLMESHFILE"]=tvmesh_file
  os.environ["SURFMESHFILE"]=tsmesh_file
  os.environ["SOURCEPOINTSFILE"]=vpmesh_file
  
  flags = " -x0 -E "
  
  if interactive:
    flags = " -e "
    
  subprocess.call(scirun_call+flags+scirun_net, shell=True)
#  os.system(scirun_call+" -E "+scirun_net )
  
  return ind_file
  
def main(argv):

  rerun = False
  interactive = False
  output_path = default_output_path
  net_path = default_net_path
  volmesh = default_volmesh
  surfmesh = default_surfmesh
  pointvol = default_sourcepoints
  tag = default_tag
  scirun = default_scirun

  opts, args = getopt.getopt(argv, "hrit:m:s:o:x:n:",
                  ["help", "rerun", "interactive", "type=", "vmesh=", "smesh=", "points=", "odir=","SCIRun=", "netpath=" ])
  
  for opt, arg in opts:
    if opt == "-h" or opt == "--help" :
      print("python RunAndCombineDipoleLF.py -r -t <source type ([pointsource], forward)> -m <tetmesh filename> -s <trisurf filename> -p <volume points filename> -o <output path> -x <scirun executable> -n <scirun network path>" )
      sys.exit()
    elif opt in ("-t", "--type"):
      tag = arg
    elif opt in ("-m", "--vmesh"):
      volmesh = arg
    elif opt in ("-s", "--smesh"):
      surfmesh = arg
    elif opt in ("-p", "--points"):
      pointvol = arg
    elif opt in ("-r", "--rerun"):
      rerun=True
    elif opt in ("-o", "--odir"):
      output_path = arg
    elif opt in ("-x", "--SCIRun"):
      scirun = arg
    elif opt in ("-n", "--netpath"):
      net_path = arg
    elif opt in ("-i", "--interactive"):
      interactive=True
      
  solution_dir = os.path.join(output_path, volmesh[:-4]+"/")
  
  if not os.path.exists(solution_dir):
    os.mkdir(solution_dir)
  
  file_root = volmesh[:-4]+param_dict[tag]["tag"]
  num_files = param_dict[tag]["fnum"]
  expectsize = param_dict[tag]["ssize"]
  scirun_net = os.path.join(net_path,param_dict[tag]["sr_net"])
  
  ind_files = sorted(glob.glob(os.path.join(solution_dir, "*"+file_root+"*index_file.txt")))
  
  for indf in ind_files:
#    print(type(indf))
#    print(indf)
    os.remove(indf)

  files, missing_files = check_files(solution_dir, file_root, num_files)
  
  if not len(files)==num_files:
    print("rerun solution for missing files")
    if len(missing_files)>0 and rerun:
      for m in missing_files:
        rerun_solution(m, file_root, scirun_net, scirun, interactive,
              os.path.join(output_path, volmesh),
              os.path.join(output_path, surfmesh),
              os.path.join(output_path, pointvol))
    return
    
  matrix = load_solutions(files, tag, expectsize)
  
#  if tag=="dipole" or tag=="pointsource":
#    LF_mat = -matrix[0].T
#  else:
#    LF_mat = matrix[0]
  
  scipy.io.savemat(os.path.join(output_path, file_root+"_LF.mat"),{"LF":matrix[0]})
  
  if len(matrix[1])==0:
    print("No problems detected")
  else:
    print(len(matrix[1]), " possible problems detected")
    if rerun:
      for f in matrix[1]:
        print("rerunning ", f)
        ind_file = rerun_solution(f, file_root, scirun_net, scirun, interactive,
              os.path.join(output_path, volmesh),
              os.path.join(output_path, surfmesh),
              os.path.join(output_path, pointvol))
              
      os.remove(ind_file)
  
  return

  
if __name__ == "__main__":
   main(sys.argv[1:])
  
  


  
  


