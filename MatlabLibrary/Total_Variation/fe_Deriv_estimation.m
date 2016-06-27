%% ESTIMATE THE FE_DERIVATIVE AND STIFFNESS MATRIX OF THE GEOMETRY


% INPUTS
ndMatH = field2.node';
elmtH = field2.cell';

conduct =eye(3);

% CODE
eleNumH = size(elemtH,1);

for i = 1:eleNumH
	
    ch = elmtH( i, 1:4);
    vx = ndMatH( ch, 1);   vy = ndMatH(ch, 2);    vz = ndMatH(ch, 3);
  
    [ LocStMat, LocDerivMat, vol ] = BuildLocStiffMat_LinearTet( 0, conduct, vx, vy, vz);
    eleVol_H( i ) = vol;
    feDerivMat( i, :, :) = LocDerivMat;
    feStiffMat( i, :, :) = LocStMat;

end


% OUTPUT
o1 = feDerivMat;
o2 = feStiffMat;

