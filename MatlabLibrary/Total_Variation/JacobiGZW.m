
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %
  %  Return the weights and zeros for Gauss-Legendre method
  %  To compute the integration of (1-x)^alpha * (1+x)^beta *u(x)
  %
  %
  %
  %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function [z,w] = JacobiGZW(degree,alpha,beta)
  
  one = 1.0; 
  two = 2.0;
  apb = alpha + beta;

  z = JacobiGZeros(degree, alpha, beta);  

  for i=1:degree,
    w(i) = JacobiPolyDerivative(degree,z(i),alpha,beta);
  end  
      
  fac  = (two^(apb + one))*gamma(alpha + ...
         degree + one)*gamma(beta + degree + one);
  fac = fac/(gamma(degree + one)*gamma(apb + degree + one));
  
  for i=1:degree,
    w(i) = fac/(w(i)*w(i)*(one-z(i)*z(i)));
  end
  w = w';
