
% Author: Rodrigo Carceroni
% Disclaimer: This code comes with no guarantee at all and its author
%   is not liable for any damage that its utilization may cause.
 
function [R,dR] = rpyMat (angs)
switch length(angs)
    case 3
% Return the 3x3 rotation matrix described by a set of Roll, Pitch and Yaw
% angles.

cosA = cos (angs(3));
sinA = sin (angs(3));
cosB = cos (angs(2));
sinB = sin (angs(2));
cosC = cos (angs(1));
sinC = sin (angs(1));

cosAsinB = cosA * sinB;
sinAsinB = sinA * sinB;

R = [ cosA*cosB  cosAsinB*sinC-sinA*cosC  cosAsinB*cosC+sinA*sinC ;
      sinA*cosB  sinAsinB*sinC+cosA*cosC  sinAsinB*cosC-cosA*sinC ;
        -sinB            cosB*sinC                 cosB*cosC         ];

if nargout > 1,
 %% also give the derivative function of R

 dR_C = [     0   cosAsinB*cosC+sinA*sinC  -cosAsinB*sinC+sinA*cosC ;
	      0  sinAsinB*cosC-cosA*sinC  -sinAsinB*sinC-cosA*cosC ;
              0            cosB*cosC                 -cosB*sinC         ];

 cosAsinB = cosA * cosB;
 sinAsinB = sinA * cosB;

 dR_B = [ -cosA*sinB  cosAsinB*sinC     cosAsinB*cosC ;
          -sinA*sinB  sinAsinB*sinC     sinAsinB*cosC ;
           -cosB        -sinB*sinC        -sinB*cosC         ];

 cosAsinB = -sinA * sinB;
 sinAsinB =  cosA * sinB;
 
 dR_A = [ -sinA*cosB  cosAsinB*sinC-cosA*cosC  cosAsinB*cosC+cosA*sinC ;
      cosA*cosB  sinAsinB*sinC-sinA*cosC  sinAsinB*cosC+sinA*sinC ;
        0            0                 0         ];
 dR_C=dR_C';
 dR_B=dR_B';
 dR_A=dR_A';
 
 dR = [dR_C(:),dR_B(:),dR_A(:)];
end
 case 1
        R=[cos(angs),sin(angs);-sin(angs),cos(angs)];
        dR=[-sin(angs);-cos(angs);cos(angs);-sin(angs)];
end   
     
end