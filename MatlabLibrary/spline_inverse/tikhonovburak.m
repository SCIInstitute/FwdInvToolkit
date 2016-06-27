function X_reg=tikhonovburak(A,L,Y,lambdas,varargin)
% Function to calculate the Greensite inverse solution
%   A - forward matrix
%   L - regularization matrix
%   Y - data
%   X_reg - inverse solution

MANUALFLAG=0;
if(numel(varargin)>0)
    if(strcmpi(varargin{1},'manual'))
        MANUALFLAG=1;
    end
end

LL=(L'*L);
AA=(A'*A);

for i=1:size(Y,2)
    
    eta=zeros(numel(lambdas),1);
    rho=zeros(numel(lambdas),1);
    for k=1:numel(lambdas)
        temp = (AA+lambdas(k)*LL)\(A'*Y(:,i));
        eta(k)=norm(L*temp);
        rho(k)=norm(Y(:,i)-A*temp);
    end
    
%     [lambda,rho_c,eta_c] = l_corner(rho,eta,lambdas);
    
    if(MANUALFLAG==1)
        figure(1), loglog(rho,eta)
        [rho_c,eta_c]=ginput(1); [val,ind]=min(sum([rho-rho_c,eta-eta_c].^2,2)); lambda=lambdas(ind);
        hold on, loglog(rho_c,eta_c,'ro'),hold off
    else
        rholog=log10(rho); etalog=log10(eta); Trho=2*numel(lambdas);
        etalog=spline(rholog,etalog,linspace(min(rholog),max(rholog),Trho));
        etalog=lowpassma(etalog,10);
        lambdaspline=spline(rholog,lambdas,linspace(min(rholog),max(rholog),Trho));
%         rholog=linspace(min(rholog),max(rholog),Trho);
        detalog=diff(etalog,2);
        signchanges=diff(sign(detalog));
        lastchangeind=find(signchanges<0,1,'last');

        lambda=lambdaspline(lastchangeind);
    end
    
%     figure(1)
%     loglog(rho,eta,'b')
%     hold on
%     loglog(10^rholog(lastchangeind),10^etalog(lastchangeind),'ro')
%     hold off
%     pause
	
	if ~numel(lambda)
		lambda = lambdas(ceil(end/2));
	end

    X_reg(:,i) = (AA+lambda*LL)\(A'*Y(:,i));
    
    fprintf('t=%i\n',i)
    
end

    
