%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Returns the zero points for Jacobi Polynomial with degree, alpha, beta
%  z is an array of dimension degree.
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function z = JacobiGZeros(degree,alpha,beta)
    if degree == 0
        z = 0;
       return
    end  
    maxit = 60;
    EPS = 1.0e-14;
    dth = pi/(2.0*degree);
 
    rlast=0.0;
    one = 1.0;
    two = 2.0;



    for k=0:degree-1,  
      r = -cos((two*k + one) * dth);
      if k
        r = 0.5*(r + rlast);
      end

      for j=1:maxit,
	
        poly = JacobiPoly(degree,r,alpha,beta);
        pder = JacobiPolyDerivative(degree,r,alpha,beta);
      
        sum = 0.0;
        for i=1:k,
          sum = sum + one/(r - z(i));
        end    
        delr = -poly / (pder - sum * poly);
        r  = r + delr;
        if abs(delr) < EPS
	  break;
        end
      end
      z(k+1)  = r;
      rlast = r;
    end
    z = z';




