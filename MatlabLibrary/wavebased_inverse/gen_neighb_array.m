function neighb_array = gen_neighb_array(pts,fac)

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

col1 = 1:size(pts,1);

neighb_array = [col1' zeros(size(pts,1),9)];
   
for i=1:size(fac,1)
    x=fac(i,1);
    y=fac(i,2);
    z=fac(i,3);
    
    % neighbors of x
    index = max(find(neighb_array(x,:)));
    y_t = 0;
    z_t = 0;
    for j = 2:index
        if(neighb_array(x,j) == y)
            y_t = 1;
        end
        if(neighb_array(x,j) == z)
            z_t = 1;
        end
    end
    if(y_t == 0)
        neighb_array(x,index+1) = y;
    end
    if(z_t == 0)
        neighb_array(x,index+2-y_t) = z;
    end
    
    %neighbors of y
    index = max(find(neighb_array(y,:)));
    x_t = 0;
    z_t = 0;
    for j = 2:index
        if(neighb_array(y,j) == x)
            x_t = 1;
        end
        if(neighb_array(y,j) == z)
            z_t = 1;
        end
    end
    if(x_t == 0)
        neighb_array(y,index+1) = x;
    end
    if(z_t == 0)
        neighb_array(y,index+2-x_t) = z;
    end
    
    %neighbors of z
    index = max(find(neighb_array(z,:)));
    x_t = 0;
    y_t = 0;
    for j = 2:index
        if(neighb_array(z,j) == x)
            x_t = 1;
        end
        if(neighb_array(z,j) == y)
            y_t = 1;
        end
    end
    if(x_t == 0)
        neighb_array(z,index+1) = x;
    end
    if(y_t == 0)
        neighb_array(z,index+2-x_t) = y;
    end
end

    
    
