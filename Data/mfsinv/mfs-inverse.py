#!/usr/bin/env python

''' Python script to implement the Method of Fundamental Solutions to solve the
	inverse problem of electrocardiology. This script is a straightforward
	implementation as outlined in:

	Y. Wang and Y. Rudy, "Application of the Method of Fundamental Solutions
	to Potential-based Inverse Electrocardiology", Ann. Biomed. Eng., 34(8)
	pp 1272-1288, 2006.

	Program written by Peter Johnston, December, 2016.

	There is one major difference between the current implementation and
	that of Wang and Rudy. Here, the ficticious points are created by
	moving the given mesh points normal to their respective surfaces.

	The program requires the following input:

	1) argv[1] - the tank file mesh (this is actually two files representing
		the closed tank surface - a points file, (x,y,z) triples, and
		a faces file (this is used to define the normal to the surface).		Note that the program assumes that the normals are outward
		pointing).

	2) argv[2] - the heart file mesh  - the comments that apply to the tank
		file mesh also apply here.

	3) argv[3] - a file of ints that maps the jacket (measurement)
		electrodes to the tank surface node numbers.
	
	4) argv[4] - a file of (x,y,z) triples that represent the locations of
		the sock electrodes. This can be the points file associated
		with argv[2].

	5) argv[5] - a file of potentials recorded at the jacket locations. 
		Bad leads can be marked with 'nan' and these potentials and
		jacket locations are adjusted accordingly.

	6) argv[6] - the value of the regularisation parameter.

	7) argv[7] - the output file containing the epicardial potentials at
		the sock electrode points.

# usage:
#
#   mfs-inverse.py tank heart electrodes sock tank_pots lambda output_file
'''

# import required routines from various packages
from sys import argv
from numpy import sqrt,pi,array,cross,zeros,einsum,arccos,sum,append,dot,diag,transpose,ones
from numpy.linalg import norm,svd
from math import isnan

# scale factor for shrinking inner surface and expanding outer surface
# scale is a dimensional unit, in this case mm.
scale_out = 40.0
scale_in = 10.0


# get information from the command line
tank_file = argv[1]
heart_file = argv[2]
electnodes = argv[3]
sock_file = argv[4]
tank_pots = argv[5]
lamb = float(argv[6])
output_file = argv[7]

# one read function
def readdata(core,dtype):
# Open data file
	datafile = open(core,'r')
# create array to hold data
	data = array([map(dtype, line.split()) for line in datafile])
# close data file
	datafile.close()
# return values
	return data

# function to write potentials
def write_epi(core,data):
    outfile = open(core+'.out','w')
    for i in range(len(data)):
        outfile.write("%f\n"%data[i])
    outfile.close()
    return

# function to normalise vectors
def normalize_v3(arr):
    ''' Normalize a numpy array of 3 component vectors shape=(n,3) '''
    lens = sqrt( arr[:,0]**2 + arr[:,1]**2 + arr[:,2]**2 )
    arr[:,0] /= lens
    arr[:,1] /= lens
    arr[:,2] /= lens                
    return arr

# function to determine length of array of vectors
def vec_len(arr):
    lens = sqrt( arr[::,0]**2 + arr[::,1]**2 + arr[::,2]**2 )
    return lens

# function to determine angle between vectors in Einstein summation notation
# (best to read numpy manual for detailed explanation)
def angle_between(a,b):
    cosang = einsum('ij,ij->i',a,b)/vec_len(a)/vec_len(b)
    return arccos(cosang)/pi

# routine to determine surface normals at mesh nodes
def norm_arr(vertices,faces):
#Create a zeroed array with the same type and shape as our vertices i.e., per
# vertex normal
    norm = zeros( vertices.shape, dtype=vertices.dtype )
    dots = zeros( faces.shape, dtype=float)
#Create an indexed view into the vertex array using the array of three indices
# for triangles
    tris = vertices[faces]
#Calculate the normal for all the triangles, by taking the cross product of
# the vectors v1-v0, and v2-v0 in each triangle             
    n = cross( tris[::,1 ] - tris[::,0]  , tris[::,2 ] - tris[::,0] )
    dots[:,0] = angle_between(tris[::,1]-tris[::,0],tris[::,2]-tris[::,0])
    dots[:,1] = angle_between(tris[::,0]-tris[::,1],tris[::,2]-tris[::,1])
    dots[:,2] = 1.0-dots[:,0]-dots[:,1]
# n is now an array of normals per triangle. The length of each normal is
# dependent the vertices, 
# we need to normalize these, so that our next step weights each normal equally.
    normalize_v3(n)
# now we have a normalized array of normals, one per triangle, i.e., per
# triangle normals.
# But instead of one per triangle (i.e., flat shading), we add to each vertex
# in that triangle, 
# the triangles' normal. Multiple triangles would then contribute to every
# vertex, so we need to normalize again afterwards.
# The cool part, we can actually add the normals through an indexed view of
# our (zeroed) per vertex normal array
#       - This bit did not seem to work, hence the for loop - PRJ.
    for i in range(len(faces)):
        norm[faces[i,0]] += n[i]*dots[i,0]
        norm[faces[i,1]] += n[i]*dots[i,1]
        norm[faces[i,2]] += n[i]*dots[i,2]
    normalize_v3(norm)
    return norm


# function to scale a surface
def scale_surf(verts,norms,scale):
    new_verts = verts+scale*norms
    return new_verts

# function to assign electrodes to node positions and normals
def assign_elec(electrodes,torso,torso_norms):
    elec_pts = zeros([len(electrodes),3])
    elec_norms = zeros([len(electrodes),3])
    for i in range(len(electrodes)):
        elec_pts[i,:] = torso[electrodes[i],:]
	elec_norms[i,:] = torso_norms[electrodes[i],:]
    return elec_pts,elec_norms

# function calculate the kernel - alternative version using vectors
def kernel(x,y):
    rad = norm(x-y)
    return 1.0/(4.0*pi*rad)

# function to calculate the normal derivative of the kernel
def normder(x,y,normal):
    rad = norm(x-y)
    numer = sum((x-y)*normal)
    return -numer/(4.0*pi*rad*rad*rad)

# routine to remove any jacket node that represents a 'bad electrode' 
# the potential for that node should be recorded as 'nan'
def cleanbspm(bspm,jacket):
    new_bspm = []
    new_jack = []
    for i in range(len(bspm)):
        if not isnan(bspm[i]):
            new_bspm.append(bspm[i])
            new_jack.append(jacket[i])
    return new_bspm,new_jack

# routine to fill coefficient matrix
# x - vector of jacket points
# y - vector of ficticious points
def get_mat(x,y,normal):
# create array of zeros of appropriate size
    mat = zeros([2*len(x),len(y)+1],dtype=float)
# set array values in the top half of the first column to 1.0
    mat[0:len(x),0] = 1.0
# loop over jacket points and ficticious points
    for i in range(len(x)):
        for j in range(len(y)):
# create upper half of coefficient matrix
            mat[i,j+1] = kernel(x[i],y[j])
# create lower half of coefficient matrix
            mat[i+len(x),j+1] = normder(x[i],y[j],normal[i])
    return mat

# routine to calculate Tikhonov inverse solution to Ax=b.
# solution assumes svd has already been done
def tik_inv(u,s,vt,b,l):
    fil_fac = s/(s*s+l*l)
    inv = dot(dot(transpose(vt),dot(diag(fil_fac),transpose(u))),b)
    return inv

# routine to calculate final epicardial potentials
#     a - vector of coefficients, 
#     x - vector required points, 
#     y - vector of ficticious pts
def epis(a,x,y):
    vec = a[0]*ones(len(x))
    mat = zeros([len(x),len(y)])
    for i in range(len(x)):
        for j in range(len(y)):
            mat[i,j] = kernel(x[i],y[j])
    vec = vec+dot(mat,a[1:])
    return vec

# start of main code

# read closed heart nodes and mesh
heart_closed_verts = readdata(heart_file+'_closed.pts',float)
heart_closed_tris = readdata(heart_file+'_closed.fac',int)
# allows for either 0 and 1 based node numbering
if heart_closed_verts.shape[0] == heart_closed_tris.max():
    heart_closed_tris = heart_closed_tris-1

# read sock nodes
sock_verts = readdata(sock_file+'.pts',float)

# read closed tank surface nodes and mesh
tank_closed_verts = readdata(tank_file+'_closed.pts',float)
tank_closed_tris = readdata(tank_file+'_closed.fac',int)
# allows for either 0 and 1 based node numbering
if tank_closed_verts.shape[0] == tank_closed_tris.max():
    tank_closed_tris = tank_closed_tris-1

# read jacket nodes
elecs = readdata(electnodes+'.nodes',int)
# allows for either 0 and 1 based node numbering
if elecs.shape[0] <= elecs.max():
    elecs = elecs-1
# flatten array to create proper vector
elecs = elecs.flatten()

# read jacket potentials
bspm = readdata(tank_pots,float)
bspm = bspm.flatten()

# clean up nan's from bspm
new_bspm,new_elecs = cleanbspm(bspm,elecs)

# get normals to vertices
inner_norms = norm_arr(heart_closed_verts,heart_closed_tris)
outer_norms = norm_arr(tank_closed_verts,tank_closed_tris)

# scale surfaces to obtain ficticious points
# note: current setup for outward pointing normals on both surfaces
new_in_verts = scale_surf(heart_closed_verts,inner_norms,-scale_in)
new_out_verts = scale_surf(tank_closed_verts,outer_norms,scale_out)

# set up vector of ficticious points
fict_pts = append(new_out_verts,new_in_verts,axis=0)

# set up rhs vector
rhs = append(new_bspm,zeros(len(new_elecs)))

# set up vector of surface points and normals at electrodes
# this assumes that the electrodes are given node numbers in the mesh.
jacket,jacket_norms = assign_elec(new_elecs,tank_closed_verts,outer_norms)

# build coefficient matrix A
amat = get_mat(jacket,fict_pts,jacket_norms)

# perform SVD of matrix A
# see numpy manual for full explanation
u,s,vt = svd(amat,full_matrices=False)

# determine u^T b
alpha = dot(transpose(u),rhs)

# calculate vector of coefficients, a, using Tikhonov
avec = tik_inv(u,s,vt,rhs,lamb)

# final epicardial potentials
epi_pot = epis(avec,sock_verts,fict_pts)

# print final epicardial potential distribution
write_epi(output_file,epi_pot)
