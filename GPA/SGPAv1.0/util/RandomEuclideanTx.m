function [Q,t]=RandomEuclideanTx(d,R,similarity)
if(nargin<3)
    similarity=0;
  if(nargin<2)
    R=1; % Generate control points in the unit hypersphere
  end
end
V=R.*randn(d,d+1);
A=V(1:d,1:d);
[Q,R]=qr(A);
alpha=det(Q);
if(similarity==0)
if(alpha<0)
    Q(:,d)=-Q(:,d);   
end
else
Q=alpha*mean(svd(R));    
end
t=V(:,d+1);

