function [x,dfdx,alpha]=steepestdescent(phi,dphidx,xinit,mingradnorm,wolfeparams)
% Algorithm: Steepest Descent
% Author: Burak Erem
% --------------------------------------
% INPUT: phi = objective function
%        dphidx = gradient of objective function
%        xinit = initialization of algorithm step sequence
%        mingradnorm = smallest 2-norm of dphidx before stopping
% --------------------------------------
% OUTPUT: x = matrix whose columns are the sequence of x points
%         dfdx = matrix whose columns are the gradients at x
%         alpha = array of the step size parameters used

% Initialization:
x=xinit;
f=phi(x);
dfdx=dphidx(x);
p=-dfdx;
k=1;
alpha(1)=1;

alphamin=1e-3;alpharate=1.5;alphamax=10;c1=1e-2;c2=0.6;
if(isfield(wolfeparams,'alphamin')),alphamin=wolfeparams.alphamin;end;
if(isfield(wolfeparams,'alphamax')),alphamax=wolfeparams.alphamax;end;
if(isfield(wolfeparams,'alpharate')),alpharate=wolfeparams.alpharate;end;
if(isfield(wolfeparams,'c1')),c1=wolfeparams.c1;end;
if(isfield(wolfeparams,'c2')),c2=wolfeparams.c2;end;

costs(k)=phi(x(:,k));
costreduction=inf;

fprintf('Init objfun=%g\n',costs(k))

% The main loop
while(costreduction>mingradnorm)

%     Line Search (satisfies strong Wolfe conditions)
    alpha(k)=linesearch(phi,dphidx,x(:,k),p(:,k),c1,c2,alphamin,alpharate,alphamax);

%     Set next sequence step x_(k+1)
    x(:,k+1)=x(:,k)+alpha(k)*p(:,k);

%     Evaluate gradient(f_(k+1))
    dfdx(:,k+1)=dphidx(x(:,k+1));

%     Set search direction p_(k+1)
    p(:,k+1)=-dfdx(:,k+1);
    
    costs(k+1)=phi(x(:,k+1));
    costreduction=costs(k)-costs(k+1);
    
%     Report progress on the current iteration
    fprintf('k=%i\t\tobjfun=%g\t\tcostreduction=%g\n',k,costs(k+1),costreduction)

%     Increment k
    k=k+1;

end
