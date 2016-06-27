  function f = Test_Right_Side_Function( problem_dimension, cond_option, conductivity, x, y , z )
  %TEST_RIGHT_SIDE,  The right side in the global Galerkin equation
  % It depends on the dimension, the way we set our test functions, and the conductivity
  
  switch problem_dimension
       %Here Test function = cos( 2*pi*x )
      case 1
          switch cond_option
              % Continuous conductivity case:  k(x) = x;
              case -1
                  f = 2*pi*sin(2*pi*x) + 4*pi*pi* (x.* cos( 2*pi*x ));
                  
              %the ECG problem we are considering about
              case 0
                  f = 0;
                  
             %conductivity = 1  on [0, 1]
             case 1
                   f = 4*pi*pi * cos( 2*pi*x );    
             
            % conductivity = 1  on [0, 0.3] & [0.6, 1],  = 3 on [ 0.3, 0.6]
            case 2
                if( conductivity == 1 )%if( (x>=0&x<=0.3) | ( x>=0.6&x<=1 ) )
                    f = 4*pi*pi * x *cos( 2*pi*x );   
                elseif( conductivity == 3 )    %elseif( x>=0.3 & x <=0.6 )
                        f = 3 * 4*pi*pi  *cos( 2*pi*x );   
                end
            
            % conductivity = 1 on [ 0, 0.4] & [0.8, 1 ],  = 0.3 on [0.4, 0.8] 
            case 3
                if( (x>=0&x<=0.4) | ( x>=0.8&x<=1 ) )
                    f = 4*pi*pi *cos( 2*pi*x );   
                elseif( x>=0.3 & x <=0.6 )
                        f = 0.3 * 4*pi*pi *cos( 2*pi*x );   
                end                
        end
        
    % 2-D problem    
    % In this case, f returned is a m*n matrix where m is the number of quadra_points in x-direction and n is the num of 
    % quadra_points in y-direction
    % Input x,y are m*n arrays each
    % Input parameter conductivity is invalid here
    case 2
        m = size(x, 1);         n = size( y, 2 ); 
        switch cond_option
          % conductivity is (x+y)*identity matrix
          case -1
            for i = 1: 1: m
              for j = 1: 1: n
                temp1 = sin( 2*pi*x(i) );
                temp2 = cos( 2*pi*x(i) );
                temp3 = sin( 2*pi*y(j) );
                temp4 = cos( 2*pi*y(j) );
                f(i,j) = 2*pi*temp1*temp4 + 2*pi*temp2*temp3 + 8*pi*pi*temp2*temp4*( x(i)+y(j) );
                          
              end
            end
          
          % case 0: the actual ECG case, no known functions
          case 0
            f = zeros(m,n);
            
          % conductivity is an identity matrix, u=cos2pix * cos2piy
          case 1
            for i = 1: 1: m
                    for j = 1:1: n
                        
                        f(i, j) = 8 * pi * pi * cos( 2*pi*x(i,j) ) * cos( 2*pi*y(i,j) );
                    end
            end
          
          % conductivity matrix is [3 0; 0 2];
          case 2
            for i = 1: 1: m
                    for j = 1:1: n
                        f(i, j) = 20 * pi * pi * cos( 2*pi*x(i) ) * cos( 2*pi*y(j) );
                    end
            end            
            
          % conductivity is identity matrix u=sin2pix*sin2piy
          case 3
            for i = 1: 1: m
                    for j = 1:1: n
                        
                        f(i, j) = 8 * pi * pi * sin( 2*pi*x(i,j) ) * sin( 2*pi*y(i,j) );
                    end
            end            
        end
        
      % 3D problem  
      case 3
          nx = size(x,1);       ny = size( y, 2);  nz = size( z, 3);
          f = zeros(nx,ny,nz);
          switch cond_option
              case 0
                  f = zeros( nx, ny, nz);
              case 1
                  for i = 1:nx
                      for j = 1:ny
                          for k = 1:nz
                            f( i,j,k ) = 12 * pi *pi * cos( 2*pi*x(i,j,k) ) * cos( 2*pi * y(i,j,k) ) * cos( 2*pi * z(i,j,k) );
                          end
                      end
                  end
                  
          end
  end
  
  return;
