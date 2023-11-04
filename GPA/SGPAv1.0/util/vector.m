function Sjvec=vector(Sj)
d=length(Sj);
Sjvec=zeros(d,d*d);
for i=1:d
Sjvec(i,:)=[zeros(1,d*(i-1)),Sj',zeros(1,d*d-d-d*(i-1))];
end
