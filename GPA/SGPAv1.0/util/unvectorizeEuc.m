function [S,A,a]=unvectorizeEuc(X,d,m,n,similarity)
if(nargin<5)
    similarity='off';
end
switch(similarity)
    case 'off'

        S=zeros(d,m);
        S(:)=X((d*d+d)*n+1:end);
        for i=[1:n]
            A{i}=zeros(d,d);
            a{i}=zeros(d,1);
            A{i}(:)=X(1+((d*d+d))*(i-1):((d*d+d))*i-d);
            a{i}(:)=X(((d*d+d))*i-d+1:((d*d+d))*i);
        end

    case 'on'
        S=zeros(d,m);
        S(:)=X((d*d+d+1)*n+1:end);
        for i=[1:n]

            A{i}=zeros(d,d);
            a{i}=zeros(d,1);
            alphai=X(1+(d*d+d+1)*(i-1));
            A{i}(:)=X(2+((d*d+d+1))*(i-1):((d*d+d+1))*i-d);
            a{i}(:)=X(((d*d+d+1))*i-d+1:((d*d+d+1))*i);
        end
end