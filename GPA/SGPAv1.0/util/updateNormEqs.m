%Update of parameters for Affine transformations with sparse
%representation.
function [deltaX,X]=updateNormEqs(X,D,lambda,V) %Levenberg-Mardquardt iteration

[d,m]=size(D{1});
n=length(D);

if nargin<4
 for i=1:n
    V{i}=ones(1,m);
 end
end

alpha=0;
for i=1:n
alpha=alpha+sum(V{i}, 2);
end

%Compute residual vectors r_{i,j}=D_{i,j}-V(S_j)A_i -a_i
r=zeros(alpha*d,1);
[S,A,a]=unvectorize(X,d,m,n);
k=1;
for i=1:n
    for j=1:m
        if(V{i}(j))
            r(d*(k-1)+1:d*k)=D{i}(:,j)-A{i}*S(:,j)-a{i};
        k=k+1;
        end
    
    end
end

%Compute U which in this case U=diag(U1,..,Um);
U=cell(1,m);
for j=1:m
U{j}=zeros(d,d);
end

for j=1:m
    for i=1:n
    if(V{i}(j))
        U{j}=U{j}+A{i}'*A{i};
    end
    end
end

%Compute Vh where Vh=diag(Vh1,...,Vhn)
Vh=cell(1,n);
for i=1:n
    Vh{i}=zeros(d*d+d,d+d*d);
    for j=1:m
        if(V{i}(j))
            Sjvec=vector(S(:,j));
            Vh{i}=Vh{i}+[Sjvec'*Sjvec,Sjvec';Sjvec,eye(d)];
        end
    end
end

%Compute W where W=[W1,W2,...,Wn]
 
for i=1:n
    W{i}=zeros(d*m,d*d+d);
    for j=1:m
        if(V{i}(j))
        W{i}(d*(j-1)+1:j*d,:)=[A{i}'*vector(S(:,j)),A{i}'];              
        end        
    end
end

%Compute epsilon_a and epsilon_b
epsilon_a=zeros(m*d,1);
for i=[1:n]
epsilon_b{i}=zeros(d*d+d,1);
end

offset=1;
k=1;
for i=[1:n]
    for j=1:m
        if(V{i}(j))
        epsilon_a(d*(j-1)+1:j*d)=epsilon_a(d*(j-1)+1:j*d)+A{i}'*r(d*(k-1)+1:k*d);
        epsilon_b{i}=epsilon_b{i}+[vector(S(:,j))'*r(d*(k-1)+1:k*d);r(d*(k-1)+1:k*d)];
        k=k+1;
        end
        
    end
end

%Compute Y1,..,Yn

H=zeros(d*m,d*m);
res=epsilon_a;
for i=[1:n]
     Vhip=Vh{i}+lambda.*eye(d*d+d,d*d+d);
    if(rcond(Vhip)<1e-12)
        disp('Rank deficient in updating normal equations. Using pseudoinverse');
        Vhipi=pinv(Vhip);
    else
        Vhipi=inv(Vhip);
    end
    
    Y{i}=W{i}*Vhipi;%inv(Vh{i}+lambda*eye(d*d+d,d+d*d));
    H=H-Y{i}*W{i}';
    res=res-Y{i}*epsilon_b{i};
end
for j=[1:m]
H(d*(j-1)+1:d*j,d*(j-1)+1:d*j)=U{j}+eye(d,d).*lambda;
end
%Compute deltaS
    
    if(rcond(H)<1e-12)
        disp('Rank deficient in updating normal equations. Using pseudoinverse');
        Hi=pinv(H);
    else
        Hi=inv(H);
    end
deltaa=Hi*res;
deltaX=zeros(size(X));

for i=[1:n]
    Vhip=Vh{i}+lambda.*eye(d*d+d,d*d+d);
    if(rcond(Vhip)<1e-12)
        disp('Rank deficient in updating normal equations. Using pseudoinverse');
        Vhipi=pinv(Vhip);
    else
        Vhipi=inv(Vhip);
    end
deltab{i}=Vhipi*(epsilon_b{i}-W{i}'*deltaa);
len=length(deltab{i});
deltaX(1+len*(i-1):i*len)=deltab{i};

end
deltaX(n*len+1:end)=deltaa;
X=X+deltaX;





