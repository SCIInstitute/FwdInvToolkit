function neighb_set=find_neighb(node_set, neighb_array)

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

%neighb_set=find_neighb(node_set, neighb_array)
    neighb_set=[];
    for i=1:length(node_set)
%        for j=1:neighb_array(node_set(i),2)
        %
        %  MY EDIT - ANDREW
        %
      
        for j = 2:max(find(neighb_array(node_set(i),:)))
            if ~ismember(neighb_array(node_set(i),j),[node_set neighb_set])
                %neighb_set=[neighb_set neighb_array(node_set(i),2+j)];
                neighb_set=[neighb_set neighb_array(node_set(i),j)];
            end
        end
    end
              
              
    
