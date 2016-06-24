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

pot_vals = [];

min_act = min(floor(activation_times))-1;
max_act = max(floor(activation_times));
length = ts.numframes;

for i = 1:238
    pi = nij_proj_pts(i,:);
    p1 = pts(nearest_pts(i,1),:);
    p2 = pts(nearest_pts(i,2),:);
    p3 = pts(nearest_pts(i,3),:);
    
    d1 = sqrt(pi*p1');
    d2 = sqrt(pi*p2');
    d3 = sqrt(pi*p3');
    dt = d1+d2+d3;
    
    a1 = floor(activation_times(nearest_pts(i,1)));
    a2 = floor(activation_times(nearest_pts(i,2)));
    a3 = floor(activation_times(nearest_pts(i,3)));
    
    ai = floor(nij_act_times(i));
    
    data1 = ts.potvals(nearest_pts(i,1),a1-min_act:a1+length-max_act);
    data2 = ts.potvals(nearest_pts(i,2),a2-min_act:a2+length-max_act);
    data3 = ts.potvals(nearest_pts(i,3),a3-min_act:a3+length-max_act);
    
    datai = [];
    for j = 1:size(data1,2)
        data_ij = data1(j)*(dt-2*d1)/dt + data2(j)*(dt-2*d2)/dt + data3(j)*(dt-2*d3)/dt;
        datai = [datai;data_ij];
    end
    
    pot_vals_i = [zeros(ai-min_act,1);datai];
    k = length-size(pot_vals_i,1);
    pot_vals_i = [pot_vals_i;zeros(k,1)];
    pot_vals = [pot_vals pot_vals_i];
end
