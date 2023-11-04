function [X]=vectorizeEuc(S,A,a,d,m,n,similarity)
if nargin<7
    similarity='off';
end

X=[];
for i=[1:n]
    switch(similarity)
        case 'off'
            X=[X;A{i}(:);a{i}];
        case 'on'
            alpha(i)=det(A{i});
            X=[X;alpha(i);A{i}(:);a{i}];
    end
end
X=[X;S(:)];
