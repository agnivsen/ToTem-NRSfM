function V=generate_missingdata(ratio,n,m,d)

%First generate how much shapes n are active per component m
nact=round(n*(1-ratio));
if(nact<2)
    disp(sprintf('Statistically not enough number of shapes for tau=%f',ratio));
    V=0;
    return;
end
singular=0;
wcond=0;
while(singular<1 || wcond==0);

for i=1:n
V{i}=zeros(m,m);
end

for j=[1:m]
    vec=randomordervec(n);
    for i=[1:nact]
    V{vec(i)}(j,j)=1;
    end
end

M=zeros(m);
for i=1:n
M=M+V{i};
end
singular=det(M);
wcond=0;
wself=1;
%Each shape must share at least d+1 points with any other
for i=[1:n]
    wself=wself&(sum(diag(V{i}).*diag(V{i}))>=(d+1));
    if(i<n)
    for k=[i+1:n]
 wcond=wcond|(sum(diag(V{i}).*diag(V{k}))>=(d+1));    
    end
    end
    end
wcond=wcond&wself;
%disp(sprintf('[%f,%f]',singular>0,wcond))
end
end

%singular=0;
%wcond=0;
%while(singular<1 || wcond==0);
%M=zeros(m);
%for i=1:n
%V{i}=diag(rand(1,m)>ratio);
%M=M+V{i};
%end
%singular=det(M);
%wcond=1;
%%Each shape must share at least d+1 points with any other
%for i=[1:n]
%%for k=[i:n]
% wcond=wcond*(sum(diag(V{i}).*diag(V{k}))>=(d+1));    
%end
%end
%wcond;
%disp(sprintf('[%f,%f]',singular>0,wcond))
%end


function vec=randomordervec(n)
ord=[1:n];
vec=zeros(1,n);
for i=[1:n] 
index=round((n-i)*rand)+1;
vec(i)=ord(index);
ord(index)=[];
end
end


