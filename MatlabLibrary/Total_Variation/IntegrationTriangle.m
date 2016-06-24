%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  This file tells how to do Integration in any 2-D triangle domain
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%function integral = IntegrationTriangle( f, vertX, vertY, z_x, weight_x, z_y, weight_y )
degree_x = 6;
degree_y = 6;
alpha1 = 0;  beta1 = 0;
alpha2 = 1;  beta2 = 0;

%%%%%%%%%%%% Build the map from the standard domain to general domain %%%%%%%%%%%%%%%%
% Here the map is from [-1,1]*[-1,1] to a general triangle domain

% The x-y coordinate for 3 vertex in each triangle mesh, vert(3) is the collapsed
vertX = [ 0 0.1 0.1 ];
vertY = [ 0 0 0.1 ];

[z_x, weight_x ] = JacobiGLZW(degree_x,alpha1,beta1);
[z_y, weight_y ] = JacobiGRZW(degree_y,alpha2,beta2);
% weight_y = weight_y/2;

num1 = size( z_x, 1 );      num2 = size( z_y, 1 );
for i = 1: 1: num1
    for j = 1: 1: num2
       x( i, j ) = vertX(1)* (1 - z_x(i))* 0.5 * (1 - z_y(j))*0.5 + vertX(2)* (1 + z_x(i))*0.5*( 1 - z_y(j) )*0.5 ...
                     + vertX(3) * ( 1 + z_y(j) )*0.5;
       y( i, j ) = vertY(1)*(1 - z_x(i)) * 0.5 * (1 - z_y(j))*0.5 + vertY(2)* (1 + z_x(i))*0.5 *(1 - z_y(j) )*0.5 ...
                     + vertY(3) * ( 1 + z_y(j) )*0.5;
    end
end
%%%%%%%%%%%% End for Mapping from standard domain to general domain%%%%%%%%%%

for i = 1: 1: num1
   for j = 1: 1: num2
      f( i, j ) = (1-10*x(i,j) )*  10*  y(i,j) ;
       %f( i, j ) = sin( x(i,j) * y(i,j) );
      %f( i, j ) = cos( x(i,j) ) * cos( y(i,j) );
      %f(i,j) = z_x(i) ^2 * z_y(j);
   end
end

% compute the partial derivative of x to the standard z_x and z_y
% x_partial_deriv( :, :,1) is the deriv to z_x,  x_partial_deriv( :, :,2) is the deriv to z_y 
 %x_partial_deriv = Differentiation_Quadrilateral( 2, x, z_x, z_y );

% compute the partial derivative of x to the standard z_x and z_y
 %y_partial_deriv = Differentiation_Quadrilateral( 2, y, z_x, z_y );

% compute the Jacobian Determinate
for i = 1: 1: num1
    for j = 1: 1: num2
%       jacob( i, j ) = x_partial_deriv( i, j, 1) * y_partial_deriv( i, j, 2 ) - x_partial_deriv( i, j, 2 ) * y_partial_deriv( i, j, 1);
      jacob(i,j) = ( vertX(1)*vertY(2)-vertX(2)*vertY(1)-vertX(1)*vertY(3)+vertX(3)*vertY(1)...
        + vertX(2)*vertY(3) - vertX(3)*vertY(2) ) / 8;
    end
end

% compute the integrand transfered from global-coordinate to standard coordinate
for i = 1: 1: num1
    for j = 1:1: num2
        integrand( i, j ) = f( i, j ) * jacob( i, j );
    end
end

integral = Integration_Quadrilateral( 2,  integrand, weight_x, weight_y);
%integral = Integration_Quadrilateral( 2,  f, weight_x, weight_y);


zMat = zeros(degree_x, degree_y );
figure(1);
mesh(z_x, z_y, zMat );
xlabel('z_x');  ylabel('z_y');

figure(2);
mesh( x, y, zMat );
xlabel('x');  ylabel('y');