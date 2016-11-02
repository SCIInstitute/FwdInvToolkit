function tau = ActGaussNewton(A,Y,L,tauinit,lambda,w,minstep)
% Implements the Gauss-Newton algorithm for solving the activation-based
% inverse problem of electrocardiography.
% => minimizes the objective function ||Y-A*X||^2+lambda*||L*X||^2 where
% X is parameterized by the C^1 polynomial approximation to a step function
% as explained in "The Depolarization Sequence of the Human Heart Surface
% Computed from Measured Body Surface Potentials" by Geertjan Huiskamp and
% Adriaan van Oosterom.
% 
% Input Variables:
% A: Forward matrix
% Y: Observations (columns index time from 1 to T=size(Y,2))
% L: Regularization matrix (typically a surface Laplacian approximation)
% lambda: Regularization parameter
% w: Width parameter in step function approximation
% tauinit: Initial phase shifts for starting the algorithm
% 
% Output Variables:
% tau: Solution phase shifts of the step functions
L=full(L);
tau=tauinit(:);
step=2*ones(size(tau));
linesteps=100; % the number steps in [0,1] to consider for the line search
alpha=linspace(0,1,linesteps);

N=size(A,2); T=size(Y,2);
u=1:T;

% if the width variable isn't an array, make it a constant array
if(length(w(:))==1),w=w*ones(N,1);end;

iter=0;

while(norm(step)>minstep)
%   Calculate the Jacobian matrix and the residuals
    J=agnjacobian(A,Y,L,tau,lambda,w);
    r=agnresidual(A,Y,L,tau,lambda,w);
    
%   Calculate the step direction
    G=J.'*J;
    while(rcond(G)<=eps)
        G=G+lambda*eye(size(G));
    end
    step=-G\J.'*r;
    
        
    
%   Perform a line search in the step direction
    err=zeros(size(alpha));
    for i=1:length(alpha)
        H=zeros(N,T);
        for n=1:N
            H(n,:)=polyactrow(u-alpha(i)*step(n)-tau(n),w(n));
        end
        err(i)=norm(Y-A*H,'fro')^2+lambda*norm(L*H,'fro')^2;
        if(i>1)
            if(err(i)>err(i-1))
                err=err(1:i);
                break;
            end
        end
    end
    
%   Update step with the result of the line search
    [val,ind]=min(err(:));
    step=alpha(ind)*step;
    
%   Update tau with the result of the step
    tau=tau+step;
    
%   Display the progress of the optimization routine
    iter=iter+1;disp(sprintf('Step:%i \t\t Size:%f',iter,norm(step)))
end

end


function r = agnresidual(A,Y,L,tau,lambda,w)

M=size(A,1); N=size(A,2); T=size(Y,2);
H=zeros(N,T);

u=1:T;
for i=1:N
    H(i,:)=polyactrow(u-tau(i),w(i));
end

E=Y-A*H;
R=sqrt(lambda)*L*H;
r=[E(:);R(:)];
end

function J = agnjacobian(A,Y,L,tau,lambda,w)

M=size(A,1); N=size(A,2); T=size(Y,2);
dE=zeros(M,T,N); dR=zeros(N,T,N);
HdH=zeros(1,T,N);

u=1:T;
for i=1:N
    HdH(1,:,i)=dpolyactrow(u-tau(i),w(i));
end

A=reshape(A,[M,1,N]);
L=reshape(L,[N,1,N]);

HdHm=repmat(HdH,[M,1,1]); % MxTxL
HdHn=repmat(HdH,[N,1,1]); % NxTxL
A=repmat(A,[1,T,1]); % MxTxL
L=repmat(L,[1,T,1]); % NxTxL

dE=A.*HdHm;
dR=-sqrt(lambda)*L.*HdHn;

J=zeros((M+N)*T,N);
for i=1:N
    tempE=dE(:,:,i); tempR=dR(:,:,i);
    J(:,i)=[tempE(:);tempR(:)];
end

end

function h = polyactrow(u,w)
u=u(:).';
h=zeros(size(u));

h(u <= -w/2)=0;
h(u >= w/2)=1;
if(w>0)
    h((-w/2 < u) & (u <= 0))=0.5*((2/w)*u((-w/2 < u) & (u <= 0))+1).^2;
    h((0 < u) & (u <= w/2))=1-0.5*((2/w)*u((0 < u) & (u <= w/2))-1).^2;
end
end

function dh = dpolyactrow(u,w)
u=u(:).';
dh=zeros(size(u));

dh(u <= -w/2)=0;
dh(u >= w/2)=1;
if(w>0)
    dh((-w/2 < u) & (u <= 0))=(4/(w^2))*u((-w/2 < u) & (u <= 0))+(2/w);
    dh((0 < u) & (u <= w/2))=-(4/(w^2))*u((0 < u) & (u <= w/2))+(2/w);
end
end

