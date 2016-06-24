%% MESSNARZ INVERSE SCRIPT FOR SCIRUN
%
%		This is a script for the matlab interface in the SCIRUN network
%		nonNegative.srn.
%
%		It implements an algorithm to solve the inverse problem in
%		electrocardiography using the method published in:
%
%			Messnarz, B., Tilg, B., Modre, R., Fischer, G., & Hanser, F. (2004). 
%			A new spatiotemporal regularization approach for reconstruction of 
%			cardiac transmembrane potential patterns. 
%			IEEE Transactions on Bio-Medical Engineering, 51(2), 273?81. doi:10.1109/TBME.2003.820394
%	
%		This method consists in a potential based inverse method that 
%		uses tikhonov reularization in a minimization algorithm that 
%		constrains the solutions to be monotonically non-decreasing.
%		
%		NETWORK INTERFACE:
%			INPUTS:
%				- i1 - double - primal/dual tradeoff parameter from the ADMM
%								algorithm.
%				- i2 - <N,T>double - ECG recordings. (N leads, T time
%									samples)
%				- i3 - <N,M>double - forward matrix.
%				- i4 - <1,3>int - lambda parameters: [minLam maxLam numLam].
%										minLam - minimum lambda = 10^minLam
%										maxLam - maximum lambda = 10^maxLam
%										numLam - number of lambdas to use
%									lambdas = 10.^linspace(minLam,maxLam,numLam);
%
%			OUTPUTS:
%				- o1 - <M,T>double - estimated solution.
%				- o2 - <M,T+1>double - dual variables of ADMM code.
%
%			DEPENDENCIES:
%				- inverse_messnarz_ADMM.m
%				- l_corner.m
%						> this piece of software is extracted from the
%						regtools software from Per Christian Hansen, IMM, July 26, 2007
%						A full version of the package is freely available in: 
%							http://www.imm.dtu.dk/~pcha/Regutools/
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

	%% DEFINES
		rho = i1;				% primal/dual tradeoff ADMM
		ECG = i2;				% BSP recordings
		A = i3;					% fwd matrix
		lambda_range = i4;		% reg param range in log10 [minLam maxLam numLam]
		margin = i5;
		
		% stopping criteria
		min_r = 1e-3;		% primal
		min_s = 1e-3;		% dual
	
	%% SET UP PROBLEM
		[T] = size(ECG,2);
		M = size(A,2);
		
		% regularization matrix
		R = eye(M);
		
		% set initial guess
		initialTMP = repmat(linspace(0,1,T),M,1);
		
		% set the regularization parameters to test
		lambda = 10.^linspace(lambda_range(1),lambda_range(2),lambda_range(3));
	
	%% APPLY INVERSE
		fprintf('Starting ADMM implementation of Messnarz inverse\n');
		tic
		for lam = 1:lambda_range(3)
			
			fprintf('Solving inverse problem for lambda %0.3f\n',log10(lambda(lam)) );
			% ADMM solver
			[TMP{lam} zk{lam}] = inverse_messnarz_ADMM(A,R,ECG,lambda(lam),initialTMP,rho,min_r,min_s,margin,false);
			logRes(lam) = (norm(ECG - A*TMP{lam},2));
			logReg(lam) = (norm(R*TMP{lam},2));
			
			% initialize next
			initialTMP = TMP{lam};
			
		end
		toc
		
	%% SELECT LAMBDA THROUGH L-CURVE
		 lambda_corner = l_corner(logRes(end:-1:1)',logReg(end:-1:1)',lambda(end:-1:1));
		
		[sink ix] =  find(lambda == lambda_corner);
		fprintf('Selected lambda: %f\n',lambda(ix));
	
	%% OUTPUT
		o1 = 100*TMP{ix}-85;	% return scaled back solution
		o2 = zk{ix};
