function [jcb, dxxi] = JcbStdTet2GeneralTet( vertX, vertY, vertZ, option)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Jacobian between std tet space  and physical tet space
% vertX, vertY: the three vertice coord of physical triangle,vert3 corresponds to collapsed
% Based on equation (4.42) in Karniadakis/Spencer book

    % x1, y1, x2, y2, x3, y3
    x1 = vertX(1);   y1 = vertY(1);     z1 = vertZ(1);
    x2 = vertX(2);   y2 = vertY(2);     z2 = vertZ(2);
    x3 = vertX(3);   y3 = vertY(3);     z3 = vertZ(3);
    x4 = vertX(4);   y4 = vertY(4);     z4 = vertZ(4);
    
    dxxi = zeros(3,3);
    switch option
        case '0to1'  % [0,1] std tet space
            dxxi(1,1) = -x1 + x2; %dx/de1
            dxxi(1,2) = -x1 + x3; %dx/de2
            dxxi(1,3) = -x1 + x4; %dx/de3

            dxxi(2,1) = -y1 + y2; %dy/de1
            dxxi(2,2) = -y1 + y3; %dy/de2
            dxxi(2,3) = -y1 + y4; %dy/de3

            dxxi(3,1) = -z1 + z2; %dz/de1
            dxxi(3,2) = -z1 + z3; %dz/de2
            dxxi(3,3) = -z1 + z4; %dz/de3
            
        case '-1to1'  % [-1,1] std tet space
            dxxi(1,1) = 0.5*(-x1 + x2); %dx/de1
            dxxi(1,2) = 0.5*(-x1 + x3); %dx/de2
            dxxi(1,3) = 0.5*(-x1 + x4); %dx/de3

            dxxi(2,1) = 0.5*(-y1 + y2); %dy/de1
            dxxi(2,2) = 0.5*(-y1 + y3); %dy/de2
            dxxi(2,3) = 0.5*(-y1 + y4); %dy/de3

            dxxi(3,1) = 0.5*(-z1 + z2); %dz/de1
            dxxi(3,2) = 0.5*(-z1 + z3); %dz/de2
            dxxi(3,3) = 0.5*(-z1 + z4); %dz/de3
    end
    jcb = det( dxxi );
return;