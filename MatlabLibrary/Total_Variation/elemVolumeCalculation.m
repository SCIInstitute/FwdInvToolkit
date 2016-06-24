function elemVolumes = elemVolumeCalculation(elem,node)
%FINDS THE VOLUMES OF THE ELEMENTS IN A MESH
%
%
%Written by: Seyhmus Guler,
%
%
%Works only for tetrahedral mesh 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%INPUTS: 
        %elem: (4 x #elements), each column containing the indices of
            %nodes of corresponding element.
        %node: (3 x #nodes), each column containing the position of
            %corresponding node. 
%OUTPUTS:
    %elemVolumes: the volumes of the elements. size: 1 x #elements    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic;    

if size(elem,1) ~= 4
    elem = elem';
end
if size(node,1) ~= 3
    node = node';
end

fprintf('Calculating the volume vector of the mesh...\n');

M = size(elem,2); %number of elements
elemVolumes = zeros(1,M);
expandedNode = reshape(node(:,elem),12,[]);

% matlabpool close force
% matlabpool open;
parfor i = 1:M
    elemVolumes(i) = abs(det([ones(1,4); reshape(expandedNode(:,i),3,4)])/6);
end
% matlabpool close;
fprintf('%s%f%s\n','The element volumes are calculated in ',toc,' seconds');



        