

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %
  % This function computes the Jocobian Polynomials of degree, alpha, beta,
  % they must be scalar
  %
  %  Input is x, x can be a either row-vector or column-vector
  % 
  %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


  function y = JacobiPoly(degree,x,alpha,beta)

  if degree == 0
    y = ones(size(x));
  elseif degree == 1
    y = 0.5*(alpha-beta+(alpha+beta+2.0)*x);
  else
    degm1 = degree-1.0; 
    tmp = 2.0*degm1+alpha+beta;
    a1= 2.0*(degm1+1)*(degm1+alpha+beta+1)*tmp;
    a2= (tmp+1)*(alpha*alpha-beta*beta);
    a3= tmp*(tmp+1.0)*(tmp+2.0);
    a4= 2.0*(degm1+alpha)*(degm1+beta)*(tmp+2.0);

    y = ((a2+a3*x).*JacobiPoly(degree-1,x,alpha,beta)- ...
	     a4*JacobiPoly(degree-2,x,alpha,beta))/a1;
  end
