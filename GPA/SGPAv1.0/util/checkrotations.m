function e=checkrotations(T)
A=T.A;
w=sign(det(A{1}));
e=1;
for i=[2:length(A)]
if(w~=sign(det(A{i})))
    e=0;
end
end

