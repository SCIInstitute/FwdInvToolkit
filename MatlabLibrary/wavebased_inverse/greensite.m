function [X_reg,lambda]=greensite(A,Y,trunc_deg)
% Function to calculate the Greensite inverse solution
%   A - forward matrix
%   Y - data
%   trunc_deg - truncation degree
%   X_reg - inverse solution

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


% set the default truncation degree to 10
if nargin < 3
    trunc_deg = 10
end

% 'whiten' the data using SVD
% use reg_tools if available
% 03/12/08 AVK - for now use built in SVD
%if exist('csvd')
%    [U,S,V]=csvd(Y'*Y);
%    [Ua,Sa,Va]=csvd(A,0);     
    
% otherwise use built in SVD
%else
    [U,S,V]=svd(Y'*Y);
    [Ua,Sa,Va]=svd(A,0);
%end

trans_Y=Y*U;
sa=diag(Sa);
Ua_trans_Y=Ua'*trans_Y;

numsa=numel(sa);

% calculate regularized solution for each frame
% until the truncation degree is reached
for i=1:trunc_deg
    % use l-curve if reg_tools is available
    if exist('l_curve')
        lambda(i)=l_curve(Ua,sa,trans_Y(:,i));
    % set default regularization parameter to .01
    else
        lambda(i) = .01;
    end
    
    X_reg(:,i)=Va(:,1:numsa)*diag(sa./(sa.^2+lambda(i)))*Ua_trans_Y(1:numsa,i);
end

% recorrelate the solution
X_reg=X_reg*U(:,1:trunc_deg)';
    
