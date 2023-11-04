function R = quat2mat(q)
% QUAT2MAT - Convert a quaternion to a 3x3 rotation matrix
%   
switch(length(q))
    case 4
a = q(1); b = q(2); c = q(3); d = q(4);
R = [a^2+b^2-c^2-d^2, 2*(b*c-a*d), 2*(b*d+a*c); ...
     2*(b*c+a*d), a^2+c^2-b^2-d^2, 2*(c*d-a*b); ...
     2*(b*d-a*c), 2*(c*d+a*b), a^2+d^2-b^2-c^2];
    case 2
    a=q(1); b=q(2);
    R=[a,b;-b,a];
end

