function h=showbody(D,V,color,linewidth,linestyle)
if(nargin<5)
    linestyle='-';
if(nargin<4)
       linewidth=1;
if(nargin<3)
    color=[1,0,0];

end
end
end
LI=[30,27;
    27,29;
    27,19;
    19,12;
    12,37;
    37,1;
    1,5;
    5,2;
    2,6;
    6,14;
    14,20;
    20,25;
    25,28;
    28,26;
    5,15;
    10,35;
    35,34;
    35,32;
    32,4;
    4,11;
    11,16;
    17,38;
    38,8;
    8,33;
    33,31;
    31,7;
    31,36;
    21,24;
    24,23;
    23,22;
    22,21;
    15,23;
    15,24;
    30,29;
    26,25;
    21,18;
    22,18;
    16,15;
    17,15;
    1,2;
    2,3;
    15,18;
    9,1;
    9,2;
    5,9;
    3,9;
    18,9;
    3,39;
    39,1;
    9,39];
k=1;
for i=[1:size(LI,1)]
if(V(LI(i,1),LI(i,1))*V(LI(i,2),LI(i,2)))
p1=D(:,LI(i,1));
p2=D(:,LI(i,2));
h(k)=line([p1(1),p2(1)],[p1(2),p2(2)],[p1(3),p2(3)]);
set(h(k),'Color',color);
set(h(k),'LineWidth',linewidth);
set(h(k),'LineStyle',linestyle);
k=k+1;
end
end

for i=[1:39]
text(D(1,i),D(2,i),D(3,i),sprintf('%d',i));
end

