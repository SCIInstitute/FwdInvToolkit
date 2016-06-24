
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %
  %
  %
  %
  %
  %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function y = JacobiPolyDerivative(degree, x, alpha, beta)

  if degree == 0
    y = zeros(size(x));
  else
    y = 0.5*(alpha+beta+degree+1)*JacobiPoly(degree-1,x,alpha+1,beta+1);
%    tmp = 2.0*degree+alpha+beta;
%    b1 = tmp*(1.0-x.*x);
%    b2 = degree*(alpha-beta-tmp*x);
%    b3 = 2.0*(degree+alpha)*(degree+beta);
%
%    y = (b2.*JacobiPoly(degree,x,alpha,beta) + ...
%        b3*JacobiPoly(degree-1,x,alpha,beta))./b1;
  end   

