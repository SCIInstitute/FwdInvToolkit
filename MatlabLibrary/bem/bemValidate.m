function bemValidate
% FUNCTION bemValidate
%
% DESCRIPTION
% This function generates some boundary element models (spherical) and computes
% both the boundary element method and the analytical solution in order to establish
% whether updates to the code still do their job
%
% INPUT -
%
% OUTPUT -
%

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

    % Model I   
    % radius 1 and 2    
    % conductivity 1

    
    model = bemGenerateSphere([2 1],[0 1 0],0.13);

    model = bemCheckModel(model);
    Transfer = bemMatrixPP(model);
    
    % Generate input and output potential distributions
    Uh = anaSolveSphere('U',model.surface{2}.pts,[0 0 0],[1 0 1],[1 2],[1 1 0]);
    Ub = anaSolveSphere('U',model.surface{1}.pts,[0 0 0],[1 0 1],[1 2],[1 1 0]);
 
    Uforward = Transfer*Uh;
    
    % Compare
    
    f1 = figure('name','original');
    bemPlotSurface(model.surface{1},Ub,'blue','colorbar');
    
    f2 = figure('name','forward BEM'); 
    bemPlotSurface(model.surface{1},Uforward,'blue','colorbar');
    
    rdm = errRDM(Ub,Uforward)
    mag = errMAG(Ub,Uforward)
