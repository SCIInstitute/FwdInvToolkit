%% HELP:
% 
% 				This function implements an algorithm that solves the inverse problem
% 				published by Messnarz et al. in TBME-IEEE:
% 					"A New Spatiotemporal Regularization Approach
% 					for Reconstruction of Cardiac Transmembrane
% 					Potential Patterns"
% 
% 				Messnarz et al. pos e least squares minimization problem
% 				defined as:
% 
% 					min_X | LX - Y |_2^2 + lambda*|RX|_2^2
% 						st. x(:,t) - x(:,t+1) <= 0
% 							x(:,1) >= minX
% 							x(:,end) <= maxX
% 
% 				where L is the forward matrix, R is the regularization matrix,
% 				X are the transmembrane potentials (TMP) and y are the
% 				ECG recordings.
% 
% 				This problem is equivalent to a quadratic problem of the form:
% 					min_x  x'Gx  +  x'c + l
% 						st. Ax >= b
% 
% 					where	G = (L'L + lambdaR'R)
% 							c = -L'Y
% 
% 				The practical implementation of this algorithm is done through
% 				the ADMM minimization of the function:
% 
% 					min_{x,z} f(x) + g(z)  + 
% 					min_{x,z} \sum_{i=1}^T \|Ax_i - y_i\|_2^2 + \sum_{j=1}^T proj( z_j ) + proj(1- z_{T+1}) 
% 
% 						st. x_{t} - x_{t-1} = z_t  ; t=[2..T]
% 							x_1 = z_1
% 							x_{T} = z_{T+1}
% 
% 						where proj is defined as: { g(z<0)=infty, g(z>=0)=0 }
% 
% 				Thus the objective the augmented lagrangian is:
% 
% 					min_{x,z,w} L(x,z,w) =
% 								\sum_{i=1}^T \|Ax_i - y_i\|_2^2 + proj( z ) + w^T( Gx(:) -b -z) + \rho/2\|Gx(:) -b - z\|_2^2
% 
% 				ADMM will solve this augmented lagrangian minimization by
% 				optimizing recursively for x, z and w (lamk in the code).
% 
% 				rho: rho determines how important is the augmented term in each
% 				iteration. It has effects on the convergence rate but and is
% 				problem dependent. Its value is updated every revisit_rho
% 				iterations if the residuals are 
%
%			INPUT:
% 					- A - <N M>double - forward operator of the linear system.
% 					- R - <H,M>double - regularization matrix.
% 					- ECG - <N,T>double - Body surface measurements.
% 					- lambda - double - regularization term.
% 					- initialx - <M,T>double - initial guess for the TMP
% 					solution.
% 					- rho - double - augmented term weight.
% 					- min_r - double - stopping criteria for the primal
% 					residual.
% 					- min_s - double - stopping criteria for the dual residual.
%					- margin - <2,1>double - minimum and maximum bounds for
%					the transmembrane potentials.
% 					- verbose - boolean - print ADMM iteratinons information.
%
%			OUTPUT:
% 					- xk - <M,T>double - estimated TMP.
% 					- zk - <M,T+1>double - estimated slack variable.
%
%			DEPENDENCES:
%
%			AUTHOR:
%					Jaume Coll-Font <jcollfont@gmail.com>
%
%
%		LICENSE:
%			The MIT License
%
%			Copyright (c) 2009 Scientific Computing and Imaging Institute,
%			University of Utah.
%
%			Permission is hereby granted, free of charge, to any person obtaining a
%			copy of this software and associated documentation files (the "Software"),
%			to deal in the Software without restriction, including without limitation
%			the rights to use, copy, modify, merge, publish, distribute, sublicense,
%			and/or sell copies of the Software, and to permit persons to whom the
%			Software is furnished to do so, subject to the following conditions:
%
%			The above copyright notice and this permission notice shall be included
%			in all copies or substantial portions of the Software.
%
%			THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
%			OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
%			FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
%			THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
%			LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
%			FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
%			DEALINGS IN THE SOFTWARE.
%
%


%% main function
function [xk,zk] = inverse_messnarz_ADMM(A,R,ECG,lambda,initialx,rho,min_r,min_s,margin,verbose)

	% DEFINE
		[N M] =size(A);
		
		revisit_rho = 100;
		
	% SET UP PROBLEM
			Q = A'*A + lambda*(R')*R;
			c = -2*A'*ECG;
			
			[U] = chol(2*Q + 2*rho*eye(M));

			
	% INITIALIZE
		xk = initialx;
		zk = [xk(:,1), (xk(:,2:end) - xk(:,1:end-1)) , xk(:,end)];
		[rk,sk] = residuals(xk,zk,zk,rho,margin);
		lamk = rho*rk;
		k = 1;
		
	% ADMM
	while true
		
		% min f(x)
			[xk] = min_L_x_overdet(U,c,xk,zk,lamk,rho);
			
		% min g(z)
			zk1 = zk;
			[zk] = min_L_z(xk,zk,lamk,rho,margin);

		% compute residuals
			[rk,sk] = residuals(xk,zk,zk1,rho,margin);
			
		% min lam
			lamk = lamk + rho*rk;
			
		% primal and dual residual norms
			nrk = norm(rk,2);
			nsk = norm(sk,2);
			
		% verbose and stopping criteria
			if verbose; fprintf('Iter: %d. Primal residual: %0.6f. Dual residual %0.6f.\n',k,nrk,nsk);end
			k = k+1;
			if ( nrk < min_r )&&( nsk < min_s )
				if verbose;fprintf('GatoDominguez!\n');end
				return;
			end
			
		% update adaptive rho
			if mod(k,revisit_rho) == 0
				[rho U] = new_rho(nrk,nsk,rho,Q,U);
			end
			
	end

end


%% min f(x) --- actual objective function (LSQ)
%		Optimize over the fitting error function. This is the Least Squares
%		problem.
%
function [xk] = min_L_x_overdet(R,c,xk,zk,lamk,rho)

	[M T] = size(xk);
	nDiv = 4;
	
	% solve for f(x_1)
	xk(:,1) = R\(R'\( -c(:,1) + rho*(zk(:,1) - zk(:,2) + xk(:,2)) - lamk(:,1) +lamk(:,2) ));
	
	% solve for f(x_T)
	xk(:,T) = R\(R'\( -c(:,T) + rho*(xk(:,T-1) + zk(:,T) + zk(:,T+1)) - lamk(:,T) - lamk(:,T+1) ));
	
	% solve for x_i i=[2:T-1] in nDiv blocks with i's sorted randomly.
	indxT = randperm(T-1);
	indxT = indxT(indxT~=1);
	for div = 1:nDiv
		indx = indxT;
		xk(:,indx) =  R\(R'\( -c(:,indx) + rho*( xk(:,indx-1) + xk(:,indx+1) + zk(:,indx) - zk(:,indx+1) ) - lamk(:,indx) +lamk(:,indx+1) ));
	end

end


%% min g(x) --- constraints
%	Optimizes over the constraint functions.
%
function [zk] = min_L_z(xk,zk,lamk,rho,margin)

	[M, T] = size(xk);
	
	% min_{z_1} g(z_1)
	zk(:,1) = xk(:,1) + 1/rho*lamk(:,1); 
	
	% min_{z_i} g(z_i) i = [2:T]
	for ii = 2:T
		zk(:,ii) = (xk(:,ii) - xk(:,ii-1)) + 1/rho*lamk(:,ii);
	end
	
	% min_{z_{T+1}} g(z_{T+1})
	zk(:,T+1) = xk(:,T) + 1/rho*lamk(:,T+1);
	
	
	% apply projections
	zk(zk(:,1)< margin(1), 1) = margin(1);
	mask = (zk(:,2:T) < 0); mask = [false(M,1) mask false(M,1)];
	zk( mask )   = 0;
	zk(zk(:,T+1)>margin(2),T+1) = margin(2);
	
end


%% compute residuals
%	Computes the new residuals (primal and dual) at each iteration.
%
function [rk,sk] = residuals(xk,zk,zk1,rho,margin)

	[M T] = size(xk);
	rk = zeros(M,T+1);
	sk = zeros(M,T+1);
	
	
	%% compute primal residuals
		rk(:,1) = xk(:,1) - zk(:,1) - margin(1);
		rk(:,2:T) = xk(:,2:end) - xk(:,1:end-1) - zk(:,2:T);
		rk(:,T+1) = xk(:,T) - zk(:,T+1) - margin(2);

	%% compute dual residuals
		dzk = rho*(zk - zk1);
		sk(:,1:T) = (dzk(:,1:T) - dzk(:,2:T+1));
		sk(:,T+1) = dzk(:,T+1);

end

%% update rho
% every revisit_rho iterations checks the difference between residuals and
% changes rho appropriately.
%
%	If		r > mu*s -> rho = tau*rho;
%	elseif	s > mu*r -> rho = 1/tau*rho;
%
function [rho,U] = new_rho(nrk,nsk,rho,Q,U)
	
	mu = 10;
	tau = 2;
	M = size(Q,1);

	if (nrk > mu*nsk)
		rho = tau*rho;
		[U] = chol(2*Q + 2*rho*eye(M));
	elseif (nsk > mu*nrk)
		rho = rho/tau;
		[U] = chol(2*Q + 2*rho*eye(M));
	end
end