%X1 and X2 must be 2xn vectors
%X3 is the closest shape to X2 and R, T are X3=R*X1+T 
function [R,T]=align_shapes(X1,X2)
n=size(X1,2);
mu1=mean(X1')';
mu2=mean(X2')';
X1n=X1;X1n(1,:)=X1n(1,:)-mu1(1);X1n(2,:)=X1n(2,:)-mu1(2);
X2n=X2;X2n(1,:)=X2n(1,:)-mu2(1);X2n(2,:)=X2n(2,:)-mu2(2);
[U1,D1,V1]=svd(X2n*X2n');
[U2,D2,V2]=svd(X1n*X1n');

R=U1*U2';
T=R*mu1-mu2;
X3=R*X1-[T(1).*ones(1,n);T(2).*ones(1,n)];