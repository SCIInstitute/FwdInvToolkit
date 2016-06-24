function points = spleval(f)
%SPLEVAL Evaluation of a spline or spline curve.
%
% points = spleval(f)
%
% Computes points on the given spline or spline curve f between
% its extreme breaks.

% Original routine fnplt by C. de Boor / latest change: Oct. 25, 1997
% Simplified by Per Christian Hansen, IMM, 04/16/98.
%
% ReguTools LICENCE
%	Copyright (c) 2008, Per Christian Hansen
%	All rights reserved.
%
%	Redistribution and use in source and binary forms, with or without
%	modification, are permitted provided that the following conditions are
%	met:
%
%		* Redistributions of source code must retain the above copyright
%		  notice, this list of conditions and the following disclaimer.
%		* Redistributions in binary form must reproduce the above copyright
%		  notice, this list of conditions and the following disclaimer in
%		  the documentation and/or other materials provided with the distribution
%
%	THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
%	AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
%	IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
%	ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
%	LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
%	CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
%	SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
%	INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
%	CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
%	ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
%	POSSIBILITY OF SUCH DAMAGE.
%
% Set default number of points.
npoints = 300;

if (f.form(1)=='B'), f = sp2pp(f); end

[breaks,coefs,l,k,d] = ppbrk(f);
x = breaks(1) + (0:npoints)*((breaks(l+1)-breaks(1))/npoints);
v=ppual(f,x);

if (d==1), points=[x;v]; else points = v; end