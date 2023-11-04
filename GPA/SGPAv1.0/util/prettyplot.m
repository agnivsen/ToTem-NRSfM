function h=prettyplot(x,y,color,marker,linestyle,method)
if(nargin<6)
    method='median';    
if(nargin<5)
method='median';    
linestyle='-';
if(nargin<4)
marker='+';
if(nargin<3)
color=[1,0,0];
end
end
end
end
if(size(y,1)>1)
stdy=std(y);

switch(method)
    case 'median'
    y=median(y);
    case 'mean'
    y=mean(y);
    case 'min'
    y=min(y);
end
end


h=plot(x(1:end),(y(1:end)),'color',color,'LineWidth',3,'marker',marker,'MarkerSize',28,'MarkerFaceColor',color,'LineStyle',linestyle);
%h=errorbar(x,y,stdy,'color',color,'LineWidth',2,'marker',marker);
