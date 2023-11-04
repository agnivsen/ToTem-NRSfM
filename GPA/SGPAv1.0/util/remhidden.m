function TRIh=remhidden(TRI,V)
TRIh=[];
for i=[1:size(TRI,1)]

if(V(TRI(i,1),TRI(i,1))*V(TRI(i,1),TRI(i,1)))
TRIh=[TRIh;TRI(i,:)];
end

end
