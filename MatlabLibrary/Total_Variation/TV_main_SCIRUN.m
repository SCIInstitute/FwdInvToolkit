%% HELP:
%
%       This script implements a Total Variation method to solve the inverse
%       problem in electrocardiography
%
%       This code is the SCIRun adaptation from the implementation of the Total Variation Method
%       published in:
%
%       Dafang Wang, Robert M. Kirby, et al. (2013).
%       "Inverse electrocardiographic source localization of ischemia: An optimization framework and finite element solution."
%       Journal of Computational Physics, Vol 250: 403-424.
%
%       This code needs of the convex optimization toolbox from CVX
%       (http://cvxr.com/cvx/)
%
%       INPUTS:
%               - field1 - struct - full bofy geometry.
%               - field2 - struct - heart geometry.
%               - field3 - struct - torso geometry.
%               - i1 - <N,1>double - ECG measured potentials.
%               - i2 - <N,1>int - indices of the geometry where ECG where
%               measured.
%               - i3 - <M,M>double - stiffness matrix on the heart.
%               - i4 - <N+M,N+M>double - stiffness matrix on the whole body.
%
%       OUTPUTS:
%               - o1 - <M,1>double - estimated inverse solution.
%
%       DEPENDECNIES:
%               - cvx - convex optimization package
%
%       AUTHOR:
%               SCIRun adaptation from:
%                       Jaume Coll-Font <jcollfont@ece.neu.edu>
%               Original algorithm and code from:
%                       Dafang Wang <dafang.wang@jhu.edu>
%
%

%% LOAD FROM SCIRUN
	% geometries
	body = field1;				% full body geometry
	heart = field2;				% heart geometry
	torso = field3;				% torso geometry
	
	% potentials
	uT = i1;					% recorded potentials on the torso
	lead_IX = i2;				% lead position indices
	
	% stiffness matrices
	stMatH_std = sparse(i3);	% heart stiffness matrix
	gStMat = sparse(i4);		% whole body stiffness matrix
	
% LOAD WITHIN MATLAB
	TV_load_SCIRUN;

	
%% DEFINE
	% startup CVX
	cvx_path = string2;
	run(strcat(cvx_path,'cvx_startup.m'));



%% CODE

cvx_solver sdpt3
epsln = 1e-6; %control the TV sharpness
beta = 1e-6; %weight for the TV term

if 1
    noiseRate = 0.001;  
    ndNumT = length(uT);
    uTRMS= norm(uT,2) / sqrt(length(uT));
    e = noiseRate* uTRMS * randn( ndNumT,1); %generate noise
    b = uT+ e;  
end
misfit = 1.5*noiseRate; %stop when the torso misfit is within this range of torso data

% % --------- set matrix -----------
% feGradMag0 = ComputeGradMagPerTet( elmtH, feDerivMat, vTMP0); %true grad mag
% feTV0 = dot( eleVol_H, feGradMag0); %true vTMP total variation

%--------- Obtain initial guess by solving a Tikhonov --
cvx_begin
    cvx_problem
    variables  vuT(ndNumT)  vu(ndNum) vH(ndNumH);
    minimize( norm(vuT-b, 2) + beta* quad_form(vH,stMatH_std) )
    subject to
        gStMat * vu == mR * vH;
        vuT == mQ*vu;
cvx_end
vTMP = vH;    clear vuT vu vH
%Otherwise, set vTMP=0:
if ~exist( 'vTMP');
    vTMP = zeros( ndNumH, 1);   
end

%% ------Start CVX iteration for Total Variation------------
iter = 0;
tol = 1e-6; %controls fix-point stop criterion
objNew = 1e10;  objOld = 1e10;
TotalTime_fixpt = [];
sID = sparse(1:ndNumH,1:ndNumH,1);

while 1 % iteration over beta
    iter = iter+1
    iterBeta = 0;
    feGradMag = ComputeGradMagPerTet( elmtH, feDerivMat, vTMP);
    objVec = [];
    vTMPold = vTMP;
    tic
    while 1 %fixed point iteration
        iterBeta = iterBeta+1
        
        dTVmat = UpdateTotalVarDeriv( elmtH, feStiffMat,  vTMP, feGradMag, epsln);
        cvx_begin
            cvx_problem
            variables  vuT(ndNumT)  vu(ndNum) vH(ndNumH);
            minimize( norm(vuT-b, 2) + beta*quad_form( vH, dTVmat+1e-13*sID) )
%alternatively:  minimize( quad_form(vuT-b, massMatT) + beta*quad_form( vH, dTVmat) )
            subject to
                gStMat * vu == mR * vH;
                vuT == mQ*vu;
        cvx_end
        vTMP = vH;
        feGradMag = ComputeGradMagPerTet( elmtH, feDerivMat, vTMP);
        feTV = dot( eleVol_H, feGradMag);
        objNew = quad_form( vuT - b, massMatT) + beta * dot( eleVol_H, feGradMag);
        objVec = [objVec; objNew]
        
        %check stop condition
        if iterBeta >= 20 fixptExitMsg='fix point reach max iteration';  break;  end
        if norm( vTMP - vTMPold )<  (1e-2 * norm(vTMP) )
            fixptExitMsg='Finish because vTMP converge';  break;
        end
        if abs(objOld - objNew) < max(8e-5, tol * norm( vTMP - vTMPold ))
            fixptExitMsg='Finish because obj improvement is small';  break;
        end
        vTMPold = vTMP;
        objOld = objNew;
        clear dTVmat vuT vu vH; cvx_clear
    end
    TotalTime_fixpt = [TotalTime_fixpt; toc];
    if norm(vuT-b) <= misfit * norm(b) %check if need to reduce beta
        break;
    elseif iter > 10
             disp('Discrepancy not attained\n');  break;
    else
        ratio = norm(vuT-b) / (misfit*norm(b))
        if ratio > 5  beta = beta / ratio
        else  beta = beta / 5; end
        clear vu vuT vH;    
    end    
end



	%% OUTPUT
		o1 = vTMP;
