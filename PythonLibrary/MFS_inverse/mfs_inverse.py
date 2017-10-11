#!/usr/bin/env python

''' Python script to implement the Method of Fundamental Solutions to solve the
	inverse problem of electrocardiology. This script is a straightforward
	implementation as outlined in:

	Y. Wang and Y. Rudy, "Application of the Method of Fundamental Solutions
	to Potential-based Inverse Electrocardiology", Ann. Biomed. Eng., 34(8)
	pp 1272-1288, 2006.

	Program written by Peter Johnston, December, 2016. adapted for FWD/INV 
    toolkit by Jess Tate, January, 2017.

	There is one major difference between the current implementation and
	that of Wang and Rudy. Here, the ficticious points are created by
	moving the given mesh points normal to their respective surfaces.


# usage:
#
#   mfs_inverse(heart,tank,sock, tank_pots, electrodes, lambda, scale_out, scale_in)


    The program requires the following input:

    1) heart - the heart triangle mesh  -   This should be a dictionary that has entries 
    for the vertices ('node') and the triangles ('face').  The data are assumed to be lists 
    of lists.  vertices are lists of x,y,z coordinates, and the faces are lists of vertices 
    for each face. Normals are assumed to be outward.

    2) tank - the tank triangle mesh of the same format as the heart surface
    
    3) sock - list of (x,y,z) triples that represent the locations of
    the sock electrodes. This can be the points file associated with the heart mesh
    
    a file of ints that maps the jacket (measurement)
    electrodes to the tank surface node numbers.

    4) tank_pots - potentials recorded at the jacket locations.
    Bad leads can be marked with 'nan' and these potentials and
    jacket locations are adjusted accordingly.

    5) electnums - list of of ints that maps the jacket (measurement)
    electrodes to the tank surface node numbers.

    6) scale_out - distance to move points of the tank surface
    
    7) scale_in - distance to move points of the heart surface
    
    8) lmbd_method - method of picking the regularization paramter
    
    9) lmbd - the value of the regularisation parameter. if lmbd_method is manual, 
    lmbd must be a single value, otherwise it should be a list of possible values,
    eg, numpy.linspace(0.0001,10,100).tolist()
    
    10) viz_regularization - (boolean) option to plot the values of the regularization values
    
    11) gammaa - addition parameter needed for rgcv method.
    
    The ouput is:
    
    epi_out - epicaridal pontentials saved as a list of lists.  
    
    
'''

# import required routines from various packages
from sys import argv
from numpy import sqrt,pi,array,cross,zeros,einsum,arccos,sum,append,dot,diag,transpose,ones,reshape,isnan,linspace,where,diff,sign,argmin,argmax
from numpy.linalg import norm,svd
#from math import isnan

def mfs_inverse(heart,tank,sock,tank_pots,electnums,scale_out,scale_in,lmbd_method,lmbda,viz_regularization,gamma):





# read closed heart nodes and mesh
    heart_closed_verts = array(heart['node'])
    heart_closed_tris = array(heart['face'])
# allows for either 0 and 1 based node numbering
    if heart_closed_verts.shape[0] == heart_closed_tris.max():
        heart_closed_tris = heart_closed_tris-1

# read sock nodes
    sock_verts = array(sock['node'])

# read closed tank surface nodes and mesh
    tank_closed_verts = array(tank['node'])
    tank_closed_tris = array(tank['face'])
# allows for either 0 and 1 based node numbering
    if tank_closed_verts.shape[0] == tank_closed_tris.max():
        tank_closed_tris = tank_closed_tris-1

# read jacket nodes
    elecs = array(electnums)
# allows for either 0 and 1 based node numbering
    if elecs.shape[0] <= elecs.max():
        elecs = elecs-1
# flatten array to create proper vector
    print(elecs.shape)
#elecs = elecs.flatten().T
    print('elecs =')
    print(elecs.shape)
    #    print(elecs)

# read jacket potentials
    bspm = array(tank_pots)
    sz_bspm = bspm.shape
    print(type(bspm))
    #bspm = bspm.flatten()
    print('bspm =')
    print(bspm.shape)
#    print(bspm)
# clean up nan's from bspm
    new_bspm,new_elecs = cleanbspm(bspm,elecs)
    print('bspm new=')
    print(bspm.shape)
#print(new_bspm)

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
    rhs = append(new_bspm,zeros((new_elecs.size,sz_bspm[1])), axis=0)
    print(rhs.shape)
#print(rhs)

# set up vector of surface points and normals at electrodes
# this assumes that the electrodes are given node numbers in the mesh.
    jacket,jacket_norms = assign_elec(new_elecs,tank_closed_verts,outer_norms)
    #print(jacket)
#print(jacket_norms)
# build coefficient matrix A
    amat = get_mat(jacket,fict_pts,jacket_norms)
    print('amat =')
    print(amat.shape)

# perform SVD of matrix A
# see numpy manual for full explanation
    u,s,vt = svd(amat,full_matrices=False)

# determine u^T b
    alpha = dot(transpose(u),rhs)
    print('alpha =')
    print(alpha.shape)
    #print(alpha)

# find regularization parameter lambda
    if lmbd_method == 'manual':
        if type(lmbda) is not float:
            raise ValueError('manual option needs a single lambda value')
        lamb = lmbda
    elif min(lmbda)<=0:
      raise ValueError('lmbda values should be greater than zero')
    elif lmbd_method == 'l_curve':
        lamb = find_lambda_lc(lmbda,u,s,vt,rhs,amat,alpha,viz_regularization)
    elif lmbd_method == 'zc':
        lamb = find_lambda_zc(lmbda,u,s,vt,rhs,amat,viz_regularization)
    elif lmbd_method == 'cresco':
        lamb = find_lambda_cresco(lmbda,s,alpha,viz_regularization)
    elif lmbd_method == 'gcv':
        lamb = find_lambda_gcv(lmbda,amat,u,s,vt,alpha,rhs,viz_regularization)
    elif lmbd_method == 'rgcv':
        lamb = find_lambda_gcv(lmbda,gamma,amat,u,s,vt,alpha,viz_regularization)
    else:
        raise ValueError('lmbd_method input must be: manual, l_curve, zc, cresco, gcv, or rgcv')

    print('running test')
    print('running lc')
    print('s shape = \f',s.shape)
    lamb = find_lambda_lc(lmbda,u,s,vt,rhs,amat,alpha,viz_regularization)
    print('running zc')
    print('s shape = \f',s.shape)
    lamb = find_lambda_zc(lmbda,u,s,vt,rhs,amat,viz_regularization)
    print('running cresco')
    print('s shape = \f',s.shape)
    lamb = find_lambda_cresco(lmbda,s,alpha,viz_regularization)
    print('running gcv')
    print('s shape = \f',s.shape)
    lamb = find_lambda_gcv(lmbda,amat,u,s,vt,alpha,rhs,viz_regularization)
    print('running rgcv')
    print('s shape = \f',s.shape)
    lamb = find_lambda_rgcv(lmbda,gamma,amat,u,s,vt,alpha,rhs,viz_regularization)


# calculate vector of coefficients, a, using Tikhonov
    avec = tik_inv(u,s,vt,rhs,lamb)
    print('avecs =')
    print(avec.shape)
#print(avec)


# final epicardial potentials
    epi_pot = epis(avec,sock_verts,fict_pts)
    print('epi_pot =')
    print(epi_pot.shape)
#print(epi_pot)

    return epi_pot


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
        elec_pts[i,:] = torso[int(electrodes[i]),:]
        elec_norms[i,:] = torso_norms[int(electrodes[i]),:]

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
    new_jack=jacket[~isnan(bspm).any(axis=1)]
    new_bspm=bspm[~isnan(bspm).any(axis=1)]
    #print(~isnan(bspm).any(axis=1))
    #print('jacket = ')
    #print(jacket)
    #print('bspm =')
    #print(bspm)

  #    new_bspm = []
  #  new_jack = []
  #  for i in range(len(bspm)):
  #      if not isnan(bspm[i]):
  #          new_bspm.append(bspm[i])
  #          new_jack.append(jacket[i])
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
    print('a=')
    print(a.shape)
    print('x=')
    print(x.shape)
    #print(len(x))
    print('y=')
    print(y.shape)
    #print(len(y))
    
    vec = a[0]*ones((len(x),1))
    #print(vec)
    mat = zeros([len(x),len(y)])
    for i in range(len(x)):
        for j in range(len(y)):
            mat[i,j] = kernel(x[i],y[j])
    vec = vec+dot(mat,a[1:])
    #    print(dot(mat,a[1:]))
#    print(vec)
    print(mat.shape)
#print(mat)
    return vec


# function to calculate CRESO function
def creso_fun(lamb,sigma,alpha):
  #print(lamb)
  #print(sigma.shape)
    #print(sigma)
    #print(alpha.shape)
    creso = sum(sigma**2*alpha**2*(sigma**2-3.0*lamb**2)/(sigma**2+lamb**2)**3)
    return -creso

# routine to calculate GCV function
def gcv_fun(lamb,amat,u,sigma,vt,alpha,rhs):
  # lambda cannot be 0
    
    avec_lss = tik_inv(u,sigma,vt,rhs,0.0)

    res = norm(dot(amat,avec_lss)-rhs)
    numer = sum(lamb**4*alpha**2/(lamb**2+sigma**2)**2)
    denom = amat.shape[0] - sum(sigma**2/(lamb**2+sigma**2))
    #    print lamb,res,numer,denom
    return (numer+res)/denom**2


# routine to RGCV function
# Note:  a call to the function GCV is required before calling this function
def rgcv_fun(lamb,gamma,amat,u,sigma,vt,alpha,rhs):
    gcv = gcv_fun(lamb,amat,u,sigma,vt,alpha,rhs)
    mu_fun = sum(sigma**4/(lamb**2+sigma**2)**2)
    return (gamma+(1.0-gamma)*mu_fun)*gcv

# routine to calculate the curvature of the L-curve
def l_curve(lamb,u,sigma,vt,rhs,amat,alpha):
    # lambda cannot be 0
    avec = tik_inv(u,sigma,vt,rhs,lamb)
    rho = norm(dot(amat,avec)-rhs)
    eta = norm(avec)
    rho2 = rho**2
    eta2 = eta**2
    #print(lamb)
    #print(min(abs(sigma)))
    #print(min(abs(alpha)))
    eta_dash = -4.0/lamb*sum(lamb**2*sigma**2/(sigma**2+lamb**2)**3*alpha**2)
    kappa = 2.0*eta2*rho2*(lamb**2*eta_dash*rho2+2.0*eta2*rho2*lamb+lamb**4*eta2*eta_dash)/\
      (eta_dash*(lamb**2*eta2**2+rho2**2)**(1.5))
    return kappa

# function to evaluate zero-crossing function
def zc_func(lamb,u,sigma,vt,rhs,amat):
    avec = tik_inv(u,sigma,vt,rhs,lamb)
    rho = norm(dot(amat,avec)-rhs)
    eta = norm(avec)
    return lamb**2*eta**2-rho**2

def find_lambda_zc(lamb,u,s,vt,rhs,amat,viz):
    zc = [zc_func(l,u,s,vt,rhs,amat) for l in lamb]
    print('zc =')
    print(zc)
    zca=array(zc)
    zc_ind=where(diff(sign(zca)))[0]
    
    
    if len(zc_ind)==0:
        raise ValueError('no zero crossing found with the lambda values provided')
    if len(zc_ind)==1:
        print('only one zero crossing found.  Consider increasing resolution of lambda values')

    zc_ind=zc_ind[0]

    print(type(zc_ind))
    print(zc_ind)
    m=(zc[zc_ind+1]-zc[zc_ind])/(lamb[zc_ind+1]-lamb[zc_ind])
#b=zc[zc_ind+1]-m*lam[zc_ind+1]
    zc_lambda=-(zc[zc_ind] - m*lamb[zc_ind])/m
    

    
    
    print('lambda value = \f',zc_lambda)

    if viz:
        print('regularization visualization not implemented yet')

    return zc_lambda

def find_lambda_cresco(lamb,s,alpha,viz):
    crsc=[creso_fun(l,s,alpha) for l in lamb]
    print('crsc =')
    print(crsc)
    ind =argmin(array(crsc))
    crsc_lambda = lamb[ind]

    print('lambda value = \f',crsc_lambda)
    
    if ind==0 or ind==(len(lamb)-1):
        print('min lambda found at the edge of the range. Consider expanding.')
  
    if viz:
        print('regularization visualization not implemented yet')
    
    return crsc_lambda

def find_lambda_lc(lamb,u,s,vt,rhs,amat,alpha,viz):
    lc= [l_curve(l,u,s,vt,rhs,amat,alpha) for l in lamb]
    print('lc =')
    print(lc)
    ind = argmax(array(lc))
    lc_lambda = lamb[ind]
    
    print('lambda value = \f',lc_lambda)
    
    if ind==0 or ind==(len(lamb)-1):
      print('min lambda found at the edge of the range. Consider expanding.')
  
    if viz:
        print('regularization visualization not implemented yet')

    return lc_lambda

def find_lambda_gcv(lamb,amat,u,s,vt,alpha,rhs,viz):
  #print('find_gcv sigma shape= \f',s.shape)
    gcv=[gcv_fun(l,amat,u,s,vt,alpha,rhs) for l in lamb]
    print('gcv= ')
    print(gcv)
    ind =argmax(array(gcv))
    gcv_lambda = lamb[ind]
  
    print('lambda value = \f',gcv_lambda)
  
    if ind==0 or ind==(len(lamb)-1):
        print('min lambda found at the edge of the range. Consider expanding.')

    if viz:
        print('regularization visualization not implemented yet')


    return gcv_lambda

def find_lambda_rgcv(lamb,gamma,amat,u,s,vt,alpha,rhs,viz):
    rgcv=[rgcv_fun(l,gamma,amat,u,s,vt,alpha,rhs) for l in lamb]
    print('rgcv= ')
    print(rgcv)
    ind = argmax(array(rgcv))
    rgcv_lambda = lamb[ind]
    
    print('lambda value = \f',rgcv_lambda)
    
    if ind==0 or ind==(len(lamb)-1):
      print('min lambda found at the edge of the range. Consider expanding.')

    if viz:
        print('regulariza}tion visualization not implemented yet')
    

    return rgcv_lambda




# start of main code


