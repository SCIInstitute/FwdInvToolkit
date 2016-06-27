% Total Var with Vogel's fixed point iteration, implemented by cvx
% minimization problem: vH = argmin norm(vuT-b,2) + beta* norm(vH,1),
%                                  subject to gStMat * vu == mR*vH
%                                                        vuT = mQ*vu
% vu: voltages throughout the volume domain
% vuT: torso surface voltage, a subset of vu. mQ is the projection operator
% vH: the heart voltage. This is what we want to solve for
% gStMat*vu== mR*vH is the finite element formulation of the original PDE
% equation

clear;
meshfile = 'tvFE_mesh1_Tr1.mat';   %load mesh and matrices

datafile = 'ts_ht1_ischm1_mesh1_Tr1_30db_FOT.mat';  %load data file

load(meshfile);  clear Ax Ay Az M1 M1A filename
load(datafile,'b','e', 'noiseRate','uT','vTMP0', 'gu','gStMat','mR','mQ','stMatH_std','vTMP');

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

% --------- set matrix -----------
feGradMag0 = ComputeGradMagPerTet( elmtH, feDerivMat, vTMP0); %true grad mag
feTV0 = dot( eleVol_H, feGradMag0); %true vTMP total variation

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

%------Start CVX iteration for Total Variation------------
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


%----- Evaluation --------------------
ratio = norm(vuT-b) / (misfit*norm(b))
check1 = norm(gStMat * vu -mR * vTMP, inf);

%Validation: compare vEXP_H with vEXP_H0,  vTMP with vTMP0, respectively
vEXP_H = vu( ndv_HtAll); %reconstructed extracell potential
vEXP_H0 = gu(ndv_HtAll); %truth from fwd computation

%RE adjusted as zero-baseline:
v1 = vTMP - mean(vTMP) + mean(vTMP0);   v0 = vTMP0;
RE = norm(v1-v0) / norm(v0);    

CCv = corrcoef( vTMP, vTMP0);
CCu = corrcoef( vEXP_H, vEXP_H0);

save scitmp.mat;