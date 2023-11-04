% Render skeleton
function h=renderskeleton(D,color,linewidth,linestyle)
if(nargin<4)
    linestyle='-';
if(nargin<3)
       linewidth=1;
if(nargin<2)
    color=[1,0,0];

end
end
end
LI=[10,1;
    1,20;
    1,2;
    1,4;
    4,5;
    2,3;
    3,13;
    5,15;
    1,11;
    11,8;
    11,6;
    6,7;
    8,9;
    7,17;
    9,19;
    ];
x=D(1,:);
y=D(2,:);
z=D(3,:);
%h=plot3(x,y,z,'ro','LineWidth',2);
%for i=[1:length(x)]
%text(x(i),y(i),z(i),sprintf('%d',i));
%end
k=1;
hold on;
for i=[1:length(x)]
        
plot3(x(i),y(i),z(i),'Marker','o','Color',color,'LineWidth',linewidth*1.5);
end
for i=[1:size(LI,1)]
p1=D(:,LI(i,1));
p2=D(:,LI(i,2));
h(k)=line([p1(1),p2(1)],[p1(2),p2(2)],[p1(3),p2(3)]);
set(h(k),'Color',color);
set(h(k),'LineWidth',linewidth);
set(h(k),'LineStyle',linestyle);
k=k+1;
end
hold off;