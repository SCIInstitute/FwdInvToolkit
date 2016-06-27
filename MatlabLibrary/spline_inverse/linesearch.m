function alphafinal = linesearch(phi,dphidx,x,p,c1,c2,alphamin,alpharate,alphamax)
% Algorithm: Line Search Algorithm (Algorithm 3.2, pg.59 of "Numerical
% Optimization" by Nocedal & Wright
% Author: Burak Erem
% ----------------------------------------
% INPUT: phi = objective function
%        dphidx = gradient of objective function
%        x = point at which we are conducting the line search
%        p = search direction
%        c1 = "sufficient decrease condition" parameter
%        c2 = "curvature condition" parameter
%        alphamin = minimum starting alpha (should be >0)
%        alpharate = multiplicative rate at which alpha grows
%                    (i.e. next=alpharate*prev)
%        alphamax = upper bound on the chosen alpha
% ----------------------------------------
% OUTPUT: alphafinal = chosen alpha

phialpha=@(a)(phi(x+a*p));
dphidalpha=@(a)(dphidx(x+a*p)'*p);

% "Exact" line search
% alphas=logspace(log10(alphamin),log10(alphamax),5);
% alphas=linspace(alphamin,alphamax,ceil(numel(x)/10));
alphas=linspace(eps,alphamax,10);
costs=zeros(size(alphas));
for i=1:numel(alphas),costs(i)=phialpha(alphas(i));end
[val,ind]=min(costs);
alphafinal=alphas(ind);
return

% Initialize the algorithm
alpha(1)=0;alpha(2)=alphamin;
i=2;

while(1)
%     Check for first terminating condition
    currphi=phialpha(alpha(i));
    if((currphi>phialpha(0)+c1*alpha(i)*dphidalpha(0))||((currphi>=phialpha(alpha(i-1)))&&(i>2)))
        alphafinal=zoom(alpha(i-1),alpha(i),phialpha,dphidalpha,c1,c2); return;
    end

%     Check for second terminating condition
    currdphidalpha=dphidalpha(alpha(i));
    if(abs(currdphidalpha)<=-c2*dphidalpha(0))
        alphafinal=alpha(i); return;
    end
    
%     Check for third terminating condition
    if(currdphidalpha>=0)
        alphafinal=zoom(alpha(i),alpha(i-1),phialpha,dphidalpha,c1,c2); return;
    end
    
%     If we haven't terminated until now, update alpha
    alpha(i+1)=alpharate*alpha(i);
    
%     Check for last terminating condition
    if(alpha(i+1)>alphamax)
        alphafinal=alphamax; return;
    end
    
    i=i+1;
end

end

% This version of zoom uses cubic interpolation
function alpha=zoom(alphalo,alphahi,phi,dphi,c1,c2)
while(1)
    % Cubic interpolation
    d1=dphi(alphalo)+dphi(alphahi)-3*(phi(alphalo)-phi(alphahi))/(alphalo-alphahi);
    d2=sqrt(max(d1^2-dphi(alphalo)*dphi(alphahi),0));
    alpha=alphalo-(alphalo-alphahi)*(dphi(alphahi)+d2-d1)/(dphi(alphahi)-dphi(alphalo)+2*d2);
    
    % Sanity checks on values of alpha
    if(alpha>alphahi),alpha=alphahi;return;end;
    if(alpha<alphalo),alpha=alphalo;return;end;
    if(abs(alphahi-alphalo)<=eps);alpha=alphalo;return;end;
    
    % Rest of zoom algorithm
    currphi=phi(alpha);
    if((currphi>phi(0)+c1*alpha*dphi(0))||(currphi>=phi(alphalo)))
        alphahi=alpha;
    else
        currdphi=dphi(alpha);
        if(abs(currdphi)<=-c2*dphi(0))
            return;
        end
        if(currdphi*(alphahi-alphalo)>=0)
            alphahi=alphalo;
        end
        alphalo=alpha;
    end
end
end