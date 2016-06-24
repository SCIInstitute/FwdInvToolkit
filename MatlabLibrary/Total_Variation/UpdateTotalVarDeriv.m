function dTVmat = UpdateTotalVarDeriv( elmtH, feStiffMat, v, feGradMag, epsln)
% update:   grad(phi_i) dot grad(phi_j)  / sqrt( grad(v)^2+epsln)
% v is not used in this func. We just need its size.

ndNumH = length(v);
dTVmat = sparse( ndNumH, ndNumH); % the deriv of TotalVar term
if ~exist( 'epsln')
    epsln = 1e-3;
end

eleNumH = size( elmtH,1);
for i = 1: eleNumH
     ch = elmtH( i, :); %size=3 for triangle, 4 for tet
     k = 1 / sqrt( (feGradMag(i))^2 + epsln );
     dTVmat( ch, ch) = dTVmat( ch, ch) + squeeze(feStiffMat(i, :, :)) * k;
end

return;