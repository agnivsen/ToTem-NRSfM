function S=RandomShape(d,m,R)
if(nargin<3)
R=1; %Default generate the reference shape in the unit hypersphere
end
S=R.*(rand(d,m));