
 

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %
  %   Return the weights and zeros for Gauss-Radau-Legendre method
  %  To compute the integration of (1-x)^alpha * (1+x)^beta *u(x)
  %
  %
  %
  %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function [z,w] = JacobiGRZW(degree,alpha,beta)

  if degree == 1
    z(1) = 0.0;
    w(1) = 2.0;
  else
    one = 1.0;
    two = 2.0;
    apb = alpha + beta; 

    z(1) = -one;
   
    z(2:degree) = JacobiGZeros(degree-1,alpha,beta+one);    
    w = JacobiPoly(degree-1,z,alpha,beta);
      
    fac  = (two^apb)*gamma(alpha + degree)*gamma(beta + degree);
    fac  = fac/(gamma(degree)*(beta + degree)*gamma(apb + degree + 1));

    w = fac*(1-z)./(w.*w);
    w(1) = w(1)*(beta+one);
 end
 z = z';
 w = w';