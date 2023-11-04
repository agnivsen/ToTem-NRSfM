function [S,A,a]=unvectorize(X,d,m,n)
S=zeros(d,m);
S(:)=X((d*d+d)*n+1:end);
for i=[1:n]
    A{i}=zeros(d,d);
    a{i}=zeros(d,1);
    A{i}(:)=X(1+(d*d+d)*(i-1):(d*d+d)*i-d);
    a{i}(:)=X((d*d+d)*i-d+1:(d*d+d)*i);
    A{i}=A{i}';
end
