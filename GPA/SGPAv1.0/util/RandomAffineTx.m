function [A,a]=RandomAffineTx(d,R)
if(nargin<2)
  R=1; % Generate control points in the unit hypersphere
end
V=R.*randn(d,d+1);
A=V(1:d,1:d);
a=V(:,d+1);

