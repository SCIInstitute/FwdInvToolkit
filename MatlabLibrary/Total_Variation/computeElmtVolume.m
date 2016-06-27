function [elmt_vol] = computeElmtVolume(geom)
%% HELP:
%		[elmt_vol] = computeElmtVolume(geom)
%			This function computes the volume of each element in the
%			provided tetrahedral mesh.
%			The underlying assumption is that the elements are tetrahedra
%	
	
	%% define
		numEl = size(geom.elem,1);
	
	%% for each element
	for ii = 1:numEl
		
		%% compute base area
			
		
		%% compute height
		
		%% compute volume
		
	end
end