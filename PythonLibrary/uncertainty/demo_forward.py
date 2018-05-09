# Shows some PCE analysis from samples:
# - generation of samples
# - construction of PCE coefficients
# - sensitivity analysis of PCE coefficients
# - comparison of PCE predictor versus truth

import itertools
import numpy as np
import os
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


import indexing, wafp, sensitivity
from test_functions import genz_oscillatory
from opolynd import opolynd_eval
from recurrence import jacobi_recurrence

# construction of our "expensive model". Here just a collection of
# not-so-complicated test functions.
Nx = 771                    # Number of "spatial" grid points (here, space is 1D)
x = np.linspace(0, 1, Nx)   # spatial grid
#print(x)
curr_dir = os.getcwd()


def expensive_model(params):
    """
      executes a scirun network to evaluated the foward problem based on changes in position.
    """
      
#    exec(open('/Users/jess/software/test_networks/uncertainty_run_torso.py').read())
    heart_geom_file = os.path.join(curr_dir,'../../Data/tikhonov_inv/25feb97_sock_closed.fld')
    torso_geom_file = os.path.join(curr_dir,'../../Data/pot_based_BEM_forward/tank.fld')
    heart_pot_file = os.path.join(curr_dir,'../../Data/tikhonov_inv/epi_pot_128_fr_170to190.mat')

    k = 10

    pots = run_network(heart_geom_file,torso_geom_file,heart_pot_file,params)
    
    return pots[:,k]

def run_network(heart_geom_file,torso_geom_file,heart_pot_file,translation_params):
  
    if not translation_params.size==3:
      raise ValueError('translation params must be length 3')
        # relative path to SCIRun executable
    scirun_call = os.path.join(curr_dir,'../../../testing_software/SCIRun/bin/SCIRun/SCIRun_test')
    scirun_net = os.path.join(curr_dir,'../../Networks/uncertainty-forward/uncertainty_forward.srn5')
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

    print( scirun_call+" -0 -x -S "+s_file_name)
    output=os.system(scirun_call+" -0 -x -S "+s_file_name)
    
    pot_solution = np.loadtxt(output_filename)
    
    os.remove(s_file_name)
    os.remove(output_filename)
    
    return pot_solution

def run_visualization(sens_output_file,sens_tot_file,sens_global_file,gt_output_file,pce_output_file):
    scirun_call = os.path.join(curr_dir,'../../../testing_software/SCIRun/bin/SCIRun/SCIRun_test')
    scirun_net= os.path.join(curr_dir,'../../Networks/uncertainty-forward/uncertainty_forward_visualization.srn5')
    torso_geom_file = os.path.join(curr_dir,'../../Data/pot_based_BEM_forward/tank.fld')

    s_file_name = torso_geom_file+'.py'
    s_file=open(s_file_name,'w+')
    s_file.write("scirun_load_network('"+scirun_net+"')\n")
    s_file.write("scirun_set_module_state('ReadField:0','Filename','"+torso_geom_file+"')\n")
    s_file.write("scirun_set_module_state('ReadMatrix:0','Filename','"+sens_output_file+"')\n")
    s_file.write("scirun_set_module_state('ReadMatrix:1','Filename','"+gt_output_file+"')\n")
    s_file.write("scirun_set_module_state('ReadMatrix:2','Filename','"+sens_tot_file+"')\n")
    s_file.write("scirun_set_module_state('ReadMatrix:3','Filename','"+sens_global_file+"')\n")
    s_file.write("scirun_set_module_state('ReadMatrix:4','Filename','"+pce_output_file+"')\n")
    
    s_file.write("scirun_execute_all()\n")
    s_file.close()

    output=os.system(scirun_call+" -S "+s_file_name+" &")

    #os.remove(s_file_name)


######## Step 1: generate samples

d = 3                       # dimension of random/parameter space
k = 2                       # polynomial degree (parameter space)
poly_space = indexing.total_degree_indices

lambdas = poly_space(d, k)
N = lambdas.shape[0]

M = N + 10  # Number of samples (expensive model runs) to perform. Must
            # be at least N, but the +10 is arbitrary.
            # Of course, large M  ===> better.

candidate_mesh_size = int(2e2)

print("Generating parameter mesh...")
z = wafp.legendre_wafp(lambdas, M=candidate_mesh_size)
z = wafp.legendre_wafp_enrichment(z, lambdas, M-N)
print(z)
# The samples are the array z, each row is a d-dimensional sample on the
# hypercube [-1,1]^d. The particular way these points are generated is
# random, so you'll get a different set of z each time this is fun, but
# the "quality" of the grid has relatively low variance from run to run.
# This takes a while, but these points may be stored for future use.



######## Step 2: run "expensive" model

u = np.zeros([Nx, M])

print("Evaluating model on mesh...")
print((z.shape))
for ind,zval in enumerate(z):
    u[:,ind] = expensive_model(zval)

print(u)

sens_output_file=os.path.join(curr_dir,'../../Data/tikhonov_inv/sensitivity_solutions.txt')
sens_tot_file=os.path.join(curr_dir,'../../Data/tikhonov_inv/total_sensitivity.txt')
sens_global_file=os.path.join(curr_dir,'../../Data/tikhonov_inv/global_sensitivity.txt')
gt_output_file=os.path.join(curr_dir,'../../Data/tikhonov_inv/ground_truth_solutions.txt')
pce_output_file=os.path.join(curr_dir,'../../Data/tikhonov_inv/pce_solutions.txt')
np.savetxt(sens_output_file,u,delimiter = ' ')

# Each column of u is a model run for a fixed parameter value

######## Step 3: compute PCE coefficients
print("Assembling PCE coefficients...")
ab = jacobi_recurrence(lambdas.max()+1, alpha=0., beta=0., probability=True)
V = opolynd_eval(z, lambdas, ab)
weights = np.sqrt(float(N)/float(M)) / np.sqrt(np.sum(V**2,axis=1))

# The PCE coefficients are computed as a weighted discrete least-squares
# estimator with the weights above.
coeffs = np.linalg.lstsq( (V.T*weights).T, (u*weights).T)[0].T

# Each row of coeffs contains PCE coefficients for a single x gridpoint.
# Each column of coeffs contains a particular PCE coefficient for all
# values of x.

######## Step 4: whatever postprocessing you want
print("Processing PCE coefficients...")
# Compute total sensitivities
total_sensitivities = sensitivity.pce_total_sensitivity(coeffs.T, lambdas, list(range(d)))
np.savetxt(sens_tot_file,total_sensitivities,delimiter = ' ')

# Compute global sensitivities
# Compute main-effect and main-interaction sensitivities
Js = [[j] for j in range(d)]
[Js.append(comb) for comb in itertools.combinations(list(range(d)), 2)]

global_sensitivities = sensitivity.pce_global_sensitivity(coeffs.T, lambdas, Js)
np.savetxt(sens_global_file,global_sensitivities,delimiter = ' ')

# Compare surrogate discrepancy at validation points
M_validation = 100
z_validation = np.random.uniform(0, 1, [100, d])

u_truth = np.zeros([Nx, M_validation])
for ind, zval in enumerate(z_validation):
    u_truth[:,ind] = expensive_model(zval)

np.savetxt(gt_output_file,u_truth,delimiter = ' ')

u_pce = np.zeros([Nx, M_validation])
V = opolynd_eval(z_validation, lambdas, ab)
u_pce = np.dot(V, coeffs.T).T

np.savetxt(pce_output_file,u_pce,delimiter = ' ')

# Compute l2- and max-errors for each grid point:
l2_error = np.sqrt(np.sum((u_truth - u_pce)**2, axis=1))/np.sqrt(N)
linf_error = np.max(np.abs(u_truth - u_pce), axis=1)

print("L2 error on validation mesh: {0:1.3e}\nMaximum error on validation mesh: {1:1.3e}".format(np.linalg.norm(l2_error)/np.sqrt(Nx), np.max(linf_error.flatten())))

run_visualization(sens_output_file,sens_tot_file,sens_global_file,gt_output_file)

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(z[:,0], z[:,1], z[:,2], 'r.')
plt.show()


