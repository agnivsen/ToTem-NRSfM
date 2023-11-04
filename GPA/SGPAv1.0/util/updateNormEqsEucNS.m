function [deltaX,X,Jt,r]=updateNormEqsEucNS(X,D,lambda,V) %Levenberg-Mardquardt iteration

[d,m]=size(D{1});
n=length(D);
if nargin<4
 for i=1:n
    V{i}=eye(m);
 end
end
alpha=0;
for i=1:n
alpha=alpha+sum(diag(V{i}));
end

%Compute residual vectors r_{i,j}=D_{i,j}-V(S_j)A_i -a_i
r=zeros(alpha*d,1);
[S,A,a]=unvectorizeEuc(X,d,m,n);
k=1;
for i=1:n
    for j=1:m
        if(V{i}(j,j))
            r(d*(k-1)+1:d*k)=D{i}(:,j)-A{i}*S(:,j)-a{i};
        k=k+1;
        end
    
    end
end

%% Compute the Jacobian
Jt=[];
dvec=(d*d-d)/2+d;%minimal parameterization
auxvec=zeros(1,dvec-d);
for i=1:n
    for j=1:m
    if(V{i}(j,j))
        
        M=[zeros(d,dvec*(i-1)),vector(A{i}*S(:,j))*diffARot(auxvec,d),eye(d),zeros(d,dvec*n-dvec*i)];
        Ms=[zeros(d,d*(j-1)),A{i},zeros(d,m*d-d*j)];
        Jt=[Jt;[M,Ms]];    
    end
    end
end

H=Jt'*Jt+lambda*eye(dvec*n+d*m);
deltaX=pinv(H)*Jt'*r;
%Additive update in S and a{i}'s
[DeltaS,DeltaA,Deltaa]=unvectorizeEucRot(deltaX,d,m,n);
S=S+DeltaS;
for i=[1:n]
a{i}=a{i}+Deltaa{i};

A{i}=DeltaA{i}*A{i};
if(det(A{i})<0)
    input('kk');
end

end

X=vectorizeEuc(S,A,a,d,m,n);
return









%Compute U which in this case U=diag(U1,..,Um);
U=cell(1,m);
for j=1:m
U{j}=zeros(d,d);
end

for j=1:m
    for i=1:n
    if(V{i}(j,j))
        U{j}=U{j}+A{i}'*A{i};
    end
    end
end

%Compute Vh where Vh=diag(Vh1,...,Vhn)
Vh=cell(1,n);
for i=1:n
    Vh{i}=zeros(d*d+d,d+d*d);
    for j=1:m
        if(V{i}(j,j))
            Sjvec=vector(S(:,j));
            Vh{i}=Vh{i}+[Sjvec'*Sjvec,Sjvec';Sjvec,eye(d)];
        end
    end
end

%Compute W where W=[W1,W2,...,Wn]
 
for i=1:n
    W{i}=zeros(d*m,d*d+d);
    for j=1:m
        if(V{i}(j,j))
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
        if(V{i}(j,j))
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
    Y{i}=W{i}*pinv(Vh{i}+lambda*eye(d*d+d,d+d*d));
    H=H-Y{i}*W{i}';
    res=res-Y{i}*epsilon_b{i};
end
for j=[1:m]
H(d*(j-1)+1:d*j,d*(j-1)+1:d*j)=U{j}+eye(d,d).*lambda;
end
%Compute deltaS
deltaa=pinv(H)*res;
deltaX=[];
for i=[1:n]
deltab{i}=pinv(Vh{i}+eye(d*d+d,d*d+d).*lambda)*(epsilon_b{i}-W{i}'*deltaa);
deltaX=[deltaX;deltab{i}];
end
deltaX=[deltaX;deltaa];
X=X+deltaX;





