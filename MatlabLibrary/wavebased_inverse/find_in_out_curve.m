function in_curve_flag=find_in_out_curve(curve_points,pts)

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

%in_curve_flag=find_in_out_curve(curve_points,pts)

in_curve_flag=ones(size(pts,1),1);
N=size(curve_points,1);
curve_points(N+1,:)=curve_points(1,:);
N=N+1;

zp=-40;  % projection plane
z0=60; % view point

curve_points_p=project_point_wise(curve_points,zp,z0);
pts_p=project_point_wise(pts,zp,z0);


for i=1:size(pts,1)
    v1= (curve_points_p-ones(N,1)*pts_p(i,:));
    nv1=sqrt(sum(v1'.*v1'));
    cos_ang=sum(v1(1:N-1,:)'.*v1(2:N,:)')./(nv1(1:N-1).*nv1(2:N));
    temp=cross(v1(1:N-1,:),v1(2:N,:));
    ang_sign=sign(temp(:,3)');
    ang(i)=sum(abs(acos(cos_ang)).*ang_sign);
end
in_curve_flag(find(abs(ang)>pi))=-1;
