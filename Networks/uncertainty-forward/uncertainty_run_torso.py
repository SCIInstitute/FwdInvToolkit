# script to run Forward simulation for uncertainty quantification.
import numpy as np
import os

def run_network(heart_geom_file,torso_geom_file,heart_pot_file,translation_params):
  
  if not translation_params.size==3:
    raise ValueError('translation params must be length 3')
  
  scirun_call = '/Users/jess/software/SCIRun/bin/SCIRun/SCIRun_test'
  scirun_net='/Users/jess/software/test_networks/uncertainty_forward.srn5'
  output_filename = heart_geom_file+'.tmp_solution.txt'
  
  t_params = translation_params.tolist()
  t_s = [ str(t) for t in t_params]
  t_string = "\\n".join(t_s)

  s_file_name = heart_geom_file+'.py'
  s_file=open(s_file_name,'w+')
  s_file.write("scirun_load_network('"+scirun_net+"')\n")
  s_file.write("scirun_set_module_state('ReadField:0','Filename','"+heart_geom_file+"')\n")
  s_file.write("scirun_set_module_state('ReadField:1','Filename','"+torso_geom_file+"')\n")
  s_file.write("scirun_set_module_state('ReadMatrix:0','Filename','"+heart_pot_file+"')\n")
  s_file.write("scirun_set_module_state('CreateMatrix:6','TextEntry','"+t_string+"')\n")
  s_file.write("scirun_set_module_state('WriteMatrix:0','Filename','"+output_filename+"')\n")
  s_file.write("scirun_set_module_state('WriteMatrix:0','FileTypeName','SimpleTextFile')\n")
  s_file.write("scirun_execute_all()\n")
  s_file.write("scirun_force_quit()\n")
  s_file.close()

  output=os.system(scirun_call+" -0 -x -S "+s_file_name)

  pot_solution = np.loadtxt(output_filename)

  os.remove(s_file_name)
  os.remove(output_filename)

  return pot_solution




