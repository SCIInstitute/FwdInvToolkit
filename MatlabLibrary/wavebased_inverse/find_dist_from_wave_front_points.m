function dist=find_dist_from_wave_front_points(interp_curve_points,pts)

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

% dist=find_dist_from_wave_front_points(curve_points,pts)

N=size(interp_curve_points,1);
M=size(pts,1);
if N==0
    dist=ones(1,M)*100;
    return;
end


heart_center=mean(pts);
heart_center=[20 20 0]; % changed for epi_endo heart
%
%  FIX!!! for different heart geometries
%
 target_pts = size(pts,1);
% v1=pts-ones(target_pts,1)*heart_center;
% v2=interp_curve_points-ones(N,1)*heart_center;

for i=1:target_pts
%     nv2=sqrt(sum(v2'.*v2'));
%     nv1=norm(v1(i,:));
%     ang=acos( (v1(i,:)*v2')./(nv1*nv2));
%     rad=.5*(nv2+nv1*ones(size(nv2)));
%     d=ang.*rad;
    temp=interp_curve_points-ones(N,1)*pts(i,:);
    d=sqrt(sum(temp'.*temp'));
        
    [dist(i), ind]=min(d);
end
