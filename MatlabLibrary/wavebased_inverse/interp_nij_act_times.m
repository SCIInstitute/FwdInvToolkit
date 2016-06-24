function nij_act_times = interp_nij_act_times(nij_proj_pts,pts,nearest_pts,activation_times);

% nij_proj_pts is projected geometry
% pts is original geometry
% nearest_pts are the triangle projected onto
% activation_times are the act times for the original data
% nij_act_times are the output activation times (linear interpolation)

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


nij_act_times = [];

for i = 1:size(nij_proj_pts,1)
    
    pi = nij_proj_pts(i,:);
    p1 = pts(nearest_pts(i,1),:);
    p2 = pts(nearest_pts(i,2),:);
    p3 = pts(nearest_pts(i,3),:);
    
    d1 = sqrt(pi*p1');
    d2 = sqrt(pi*p2');
    d3 = sqrt(pi*p3');
    dt = d1+d2+d3;
    
    a1 = activation_times(nearest_pts(i,1));
    a2 = activation_times(nearest_pts(i,2));
    a3 = activation_times(nearest_pts(i,3));
   
    %v = [a1;a2;a3];
    
    %vi = interp3(x,y,z,v,xi,yi,zi);
    
    vi = a1*(dt-2*d1)/dt + a2*(dt-2*d2)/dt + a3*(dt-2*d3)/dt;
    
    nij_act_times = [nij_act_times;vi];
end

    
    
