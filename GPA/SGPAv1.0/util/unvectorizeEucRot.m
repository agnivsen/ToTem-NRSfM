function [S,A,a]=unvectorizeEucRot(X,d,m,n,similarity)
if(nargin<5)
    similarity='off';
end

switch(similarity)

    case 'off'
        S=zeros(d,m);
        S(:)=X(((d*(d+1))/2)*n+1:end);
        for i=[1:n]
            A{i}=zeros(d,d);
            a{i}=zeros(d,1);
            A{i}(:)=rpyMat(X(1+((d*d+d)/2)*(i-1):((d*d+d)/2)*i-d));
            a{i}(:)=X(((d*d+d)/2)*i-d+1:((d*d+d)/2)*i);
        end
    case 'on'

        S=zeros(d,m);
        S(:)=X(((d*(d+1))/2+1)*n+1:end);
        for i=[1:n]
            A{i}=zeros(d,d);
            a{i}=zeros(d,1);
            alphai=X(1+((d*d+d)/2+1)*(i-1))+1;
            A{i}(:)=(alphai)*rpyMat(X(2+((d*d+d)/2+1)*(i-1):((d*d+d)/2+1)*i-d));
            a{i}(:)=X(((d*d+d)/2+1)*i-d+1:((d*d+d)/2+1)*i);
        end

end