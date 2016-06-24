function outtensor = tensorRightMatrixMultiply(intensor,rightIndex,rightMatrix)
% Author: Burak Erem

blockmat=num2cell(intensor,[1,rightIndex+1]);
matsize=size(blockmat{1});
cellsize=size(blockmat);
blockvecmat=squeeze(blockmat); blockvecmat=blockvecmat(:);
vecmat=squeeze(cat(1,blockvecmat{:}));

multvecmat=vecmat*rightMatrix;

REblockvecmat=mat2cell(multvecmat,matsize(1)*ones(1,size(multvecmat,1)/matsize(1)),size(multvecmat,2));
multcellsize=cellsize(setdiff(1:numel(cellsize),[1,rightIndex+1]));
REblockvecmat=reshape(REblockvecmat,[1,1,multcellsize]);
outtensor=cell2mat(REblockvecmat);
multinds=1:numel(cellsize); multinds(2:rightIndex)=3:rightIndex+1; multinds(rightIndex+1)=2;
outtensor=permute(outtensor,multinds);

