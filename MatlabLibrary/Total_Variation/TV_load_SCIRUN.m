%% HELP:
%	This scrpt loads and/or precomputes the input variales into 
%	The Total Variation Code for SCIRUN.
%
		
	%% DEFINES
		% names for input files
		string_heart_ele_vol = 'eleVol_H';
		string_heart_der = 'fe_derivMTRX';
		string_stiffness_mtrx = 'stiffnessMTRX';
		string_mass_mtrx = 'massMTRX';
		string_loc_mass_mtrx = 'mR';
		string_fe_stif_mtrx = 'feStiffMat';

		% load geometry fields into corresponding variables
		elmtH = heart.cell';	ndMatH = heart.node';		% heart
		elmtT = torso.cell';	ndMatT = torso.node';		% torso
		Elmt = body.cell';	NdMat = body.node';			% whole body
		
		
		% define sizes
		ndNum = size(NdMat,1);
		ndNumH = size(ndMatH,1);
		eleNumH = size(elmtH,1);
		eleNum = size(Elmt,1);
		
	%% HEART ELEMENTS VOLUME
			fprintf('Computing heart elements volume...');
			eleVol_H = zeros(eleNumH,1);
			for ii = 1:eleNumH
				eleVol_H(ii) = elemVolumeCalculation([1 2 3 4], ndMatH(elmtH(ii,:),:));
			end
	
	%% DERIVATIVE ESTIMATOR ON HEART
			fprintf('Computing FE derivative estimation on the heart...');
			feDerivMat = zeros(eleNumH,3,4);
			for ii = 1:eleNumH
				
				ch = elmtH(ii,:);
				
				[ LocStMat, LocDerivMat ] = BuildLocStiffMat_LinearTet( 0, eye(3,3), ndMatH( ch, 1), ndMatH(ch, 2), ndMatH(ch, 3));
				feDerivMat( ii, :, :) = LocDerivMat;
				
			end

	%% CREATE ELECTRODE SELECTION MATRIX
		nLeads = numel(lead_IX);
		mQ = spalloc(nLeads,ndNum,nLeads);
		for ii = 1:nLeads
			mQ(ii,lead_IX(ii)) = 1;
		end
		
	%% MASS MATRIX
			fprintf('Calculaing mass matrix...');
			
			[z_x, w_x] = JacobiGLZW( 4, 0, 0);
			[z_y, w_y] = JacobiGRZW( 4, 1, 0 );     w_y = w_y / 2;
			[z_z, w_z] =  JacobiGRZW( 4, 2, 0 );     w_z = w_z / 4;
			zpTet = [ z_x, z_y, z_z];          wsTet = [w_x, w_y, w_z];
			qdTet = Transform2StdTetSpace( z_x, z_y, z_z, '0to1');
			coefmatBaseTet = LocalCanonicalTetBasis('Base', 1);
			phiTet = EvalBasisTet(1, coefmatBaseTet, qdTet);

			gMassMat = sparse( ndNum, ndNum);
			for i = 1:eleNum
				cn = Elmt( i, :);
				vx = NdMat(cn, 1);  vy = NdMat(cn,2);       vz = NdMat(cn, 3);
				 [ locMassMat, jcbInt ] = BuildLocMassMat_LinearTet( 0,vx, vy, vz, wsTet, phiTet );
				 gMassMat( cn, cn) = gMassMat(cn, cn) + locMassMat;
			end
			
			massMatT = full(gMassMat( lead_IX, lead_IX));
		
	%% FE STIFFNESS MATRIX
			fprintf('Calculaing FE stiffness matrix...');
			conduct = eye(3,3);
			for i = 1:eleNumH
				ch = elmtH( i, 1:4);
				vx = ndMatH( ch, 1);   vy = ndMatH(ch, 2);    vz = ndMatH(ch, 3);

				[ LocStMat, LocDerivMat, vol ] = BuildLocStiffMat_LinearTet( 0, conduct, vx, vy, vz);
				eleVol_H( i ) = vol;
				feDerivMat( i, :, :) = LocDerivMat;
				feStiffMat( i, :, :) = LocStMat;
			%     stH( ch, ch) = stH(ch, ch) + LocStMat;
			end
			
	%% heart local mass matrix
			fprintf('Calculating local mass matrix...');
			[D] = PairwiseDistance(NdMat',ndMatH');
			[sink heart_ix] = min(D,[],1);
			mR = sparse( ndNum, ndNumH);
			mR( heart_ix, :) = (-1) *stMatH_std;
