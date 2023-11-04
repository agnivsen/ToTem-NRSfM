% match_shapes_similarity

function [R,T] = match_shapes_similarity(X1,X2)

%X : [x y]
m = size(X1,1);

% Centrer donn�es
mX1 = mean(X1,1);
mX2 = mean(X2,1);

X1c = X1 - ones(m,1)*mX1;
X2c = X2 - ones(m,1)*mX2;

% 
normX1 =  sqrt(X1c(:,1)' * X1c(:,1) + X1c(:,2)' * X1c(:,2));

a = ( X1c(:,1)' * X2c(:,1) + X1c(:,2)' * X2c(:,2) ) / normX1^2;

b = (X1c(:,1)' * X2c(:,2) -  X2c(:,1)' * X1c(:,2)) / normX1^2;

theta = atan2(b,a);

R=[cos(theta),-sin(theta);sin(theta),cos(theta)];
T=R*mX1'-mX2';

%s = sqrt(a^2+b^2);




% Calculs faits pour les donn�es centr�es. Pour recomposer on tiendra
% compte des mXi retourn�s.