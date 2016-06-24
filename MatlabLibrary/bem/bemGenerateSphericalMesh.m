function surface = bemGenerateSphericalMesh(D);

% FUNCTION surface = bemGenerateSphericalMesh(D);
%
% DESCRIPTION
% This function generates a sphere and creates a triangulted surface
% D is a measure for the resolution of the sphere
% The sphere being generated is a unity sphere

%
%  For more information, please see: http://software.sci.utah.edu
% 
%  The MIT License
% 
%  Copyright (c) 2009 Scientific Computing and Imaging Institute,
%  University of Utah.
% 
%  
%  Permission is hereby granted, free of charge, to any person obtaining a
%  copy of this software and associated documentation files (the "Software"),
%  to deal in the Software without restriction, including without limitation
%  the rights to use, copy, modify, merge, publish, distribute, sublicense,
%  and/or sell copies of the Software, and to permit persons to whom the
%  Software is furnished to do so, subject to the following conditions:
% 
%  The above copyright notice and this permission notice shall be included
%  in all copies or substantial portions of the Software.
% 
%  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
%  OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
%  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
%  THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
%  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
%  FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
%  DEALINGS IN THE SOFTWARE.
%


% Create the sphere itself

r = ceil(2/D);
dr = 2/r;

R = -1:dr:1;
Z = sin(R*pi/2);

k = 0;

for p = 1:size(R,2),

Rxy = (sqrt(1-Z(p)*Z(p)));
no = ceil(2*pi*Rxy/D); 
if no == 0, no = 1; end;
phi = 0:2*pi/no:(1-1/no)*2*pi ;
phi = phi + mod(p,2)*pi/no;

clear Slice;

for q = 1:size(phi,2),

x = Rxy*cos(phi(q));
y = Rxy*sin(phi(q));
z = Z(p);
k = k + 1;
Pos(1:3,k) = [x y z]';
Slice(q) = k;
end;
Slices{p} = Slice;
end;


% Here the triangulation starts

N = size(R,2);

Tri = zeros(3,2);
k = 1;
  
for q = 1:(N-1),

H1 = Slices{q};
H2 = Slices{q+1};

% I1 = 1;

% P = Pos(:,H1(I1))*ones(1,size(H2,2))-Pos(:,H2);
 
% L = sqrt(sum(P.*P));
% I2 = find(L  == min(L));

% if I2 > 1,
% H2 = [H2(I2:size(H2,2)) H2(1:(I2-1))]; 
% end;


I1 = 1; I2 = 1;

H1 = [H1 H1(1)];
H2 = [H2 H2(1)];
N1 = size(H1,2);
N2 = size(H2,2);

if N1 == 2, H1 = Slices{q}; N1 = 1; end;

if N2 == 2, H2 = Slices{q+1}; N2 = 1; end;

while ~((I1 == N1) & (I2 == N2)),

if (I1 < N1) & (I2 < N2)
L1 = sqrt(sum((Pos(:,H1(I1+1))-Pos(:,H2(I2))).*(Pos(:,H1(I1+1))-Pos(:,H2(I2)))));
L2 = sqrt(sum((Pos(:,H1(I1))-Pos(:,H2(I2+1))).*(Pos(:,H1(I1))-Pos(:,H2(I2+1)))));
if L1 < L2,
 Tri(:,k) = [H1(I1) H1(I1+1) H2(I2)]';
 k = k + 1;
 I1 = I1 + 1;
else
 Tri(:,k) = [H2(I2) H2(I2+1) H1(I1)]';
 k = k + 1;
 I2 = I2 + 1;
end;
end;

if (I2 == N2)
 Tri(:,k) = [H1(I1) H1(I1+1) H2(I2)]';
 k = k + 1;
 I1 = I1 + 1;
elseif (I1 == N1)
 Tri(:,k) = [H2(I2) H2(I2+1) H1(I1)]';
 k = k + 1;
 I2 = I2 + 1;
end;
  
end;
end;

R = 1/3*(Pos(:,Tri(1,:))+Pos(:,Tri(2,:))+Pos(:,Tri(3,:)));
y1 = Pos(:,Tri(1,:));
y2 = Pos(:,Tri(2,:));
y3 = Pos(:,Tri(3,:));

n = cross(y2-y1,y3-y2);
nR = sum(n.*R);
H = find(nR < 0);
Temp = Tri(2,H);
Tri(2,H) = Tri(3,H);
Tri(3,H) = Temp;

% build resulting cell structure


fprintf(1,'Generated a shell of %d nodes and %d triagles\n',size(Pos,2),size(Tri,2));

surface.node = Pos;
surface.face = Tri;

% end
