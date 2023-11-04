% Experiment with the 2D Face Dataset


% load ../datasets/2Dfacedatabase/AllShapesManu
% n=size(AllShapes,2);
% m=size(AllShapes,1);
% d=2;
% k=1;
% for i=[1:15:n]
% x=AllShapes(1:m/2,i);
% y=AllShapes(m/2+1:m,i);
% D{k}=[x';y'];
% alpha=pi.*randn;
% Rp=[cos(alpha),sin(alpha);-sin(alpha),cos(alpha)];
% %D{k}=Rp*D{k}+10.*randn(size(D{k}));
% V{k}=eye(m/2,m/2);
% frame{k}=sprintf('deform%.4d.jpg',i);
% k=k+1;
% end
% n=size(D,2);
% V=generate_missingdata(0.1,n,m/2,d);

% Load shapes
%load ../../datasets/2DfacedatabaseUAH/Shapes
load ../../datasets/2Dfacedatabase/ShapeData
n=size(D,2);
m=size(D{1},2);
d=size(D{1},1);
methods={'AFF-ALL','EUC-ALL','EUC-ALT','SIM-ALL','SIM-ALT'}
nmethods=length(methods);
Result=cell(1,nmethods);
%V=generate_missingdata(0.1,n,m,d);
% GPA with all methods
for i=1:nmethods
options.method=methods{i};
options.verbose='on';
Result{i}=gpa(D,V,options);
end

%% Plot images with overlaid shapes.
paleta=hsv(nmethods);
harray=zeros(1,nmethods+1);
fontsize=30;
for i=[1:size(D,2)]
    %I=imread(sprintf('../datasets/2DfacedatabaseUAH/%s',frame{i}));
    %I=imread(sprintf('../datasets/2Dfacedatabase/%s',frame{i}));
    I=imread(framelist{i});
    h=figure(i);
    imshow(I);hold on;
    TRI=delaunay(D{i}(1,:),D{i}(2,:));
        TRIh=remhidden(TRI,V{i});        
    ht=triplot(TRIh,D{i}(1,:),D{i}(2,:));
    hold on
    hf=get(h,'CurrentAxes');
    set(hf,'FontSize',round(fontsize/2));
    set(hf,'YDir','Reverse');
    set(ht,'Color',[0.9,0.9,0.9]);
    set(ht,'LineWidth',4,'LineStyle','-');
        ht=triplot(TRIh,D{i}(1,:),D{i}(2,:));
    hold on
    hf=get(h,'CurrentAxes');
    set(hf,'YDir','Reverse');
    set(ht,'Color',[0,0,0]);
    set(ht,'LineWidth',2,'LineStyle','-');    
    harray(1)=ht(1);
namearray{1}=sprintf('Shape %d',i);
hold off
axis([100,500,80,450]);
%print(h,'-depsc',sprintf('figs/exp2DFaceDatasetShape%d',i));
saveas(h,sprintf('figs/exp2DFaceDataset/exp2DFaceDatasetShape%d',i),'fig');
print(h,'-depsc',sprintf('figs/exp2DFaceDataset/exp2DFaceDatasetShape%d',i));
end

%% Plot images with overlaid shapes.
paleta=0.8.*[1,0,0;0,1,0;0,0,1;1,1,0;0,1,1;0,0,0;0,1,0];
linestyle={'-','-.','--','-','-.',':','--',':'};
harray=zeros(1,nmethods+1);
fontsize=30;
for i=[1:size(D,2)]
    %I=imread(sprintf('../datasets/2DfacedatabaseUAH/%s',frame{i}));
    %I=imread(sprintf('../datasets/2Dfacedatabase/%s',frame{i}));
    %I=imread(framelist{i});
    h=figure(i);
   % imshow(I);hold on;
    TRI=delaunay(D{i}(1,:),D{i}(2,:));
        TRIh=remhidden(TRI,V{i});        
    ht=triplot(TRIh,D{i}(1,:),D{i}(2,:));
    hold on
    hf=get(h,'CurrentAxes');
    set(hf,'YDir','Reverse');
    set(ht,'Color',[0,0,0]);
    set(ht,'LineWidth',2,'LineStyle','--');
harray(1)=ht(1);
namearray{1}=sprintf('Shape %d',i);
for j=1:nmethods
    S=Result{j}.S;
    A=Result{j}.A;
    a=Result{j}.a;
    Di=A{i}*S+a{i}*ones(1,m);
      ht=triplot(TRI,Di(1,:),Di(2,:));
    hold on
    harray(j+1)=ht(1);
    namearray{j+1}=methods{j};
    hf=get(h,'CurrentAxes');
        set(hf,'FontSize',round(fontsize/2));
    set(hf,'YDir','Reverse');
    set(ht,'Color',paleta(j,:));
    set(ht,'LineWidth',2,'LineStyle',linestyle{j});
end
if(i==1)
hl=legend(harray,namearray);
   set(hl,'FontSize',fontsize,'Location','NorthEast');
end
   hold off
%print(h,'-depsc',sprintf('figs/exp2DFaceDatasetShape%d',i));
saveas(h,sprintf('figs/exp2DFaceDataset/exp2DFaceDatasetErrors%d',i),'fig');
print(h,'-depsc',sprintf('figs/exp2DFaceDataset/exp2DFaceDatasetErrors%d',i));
end


%% Plot reference shapes for all methods.

for j=[1:nmethods]
h=figure(j+n);
hold on;
harray=zeros(1,2);
[error,T2]=referencespaceerror(Result{j},D,V); 
S=T2.S;
    A=T2.A;
    a=T2.a;
    for i=1:n
    Dip=inv(A{i})*D{i}-inv(A{i})*a{i}*ones(1,m);
    ht=triplot(TRI,Dip(1,:),Dip(2,:));
    harray(2)=ht(1);
    hf=get(h,'CurrentAxes');
    set(hf,'YDir','Reverse');
    set(ht,'Color',[0,0,1]);
    set(ht,'LineWidth',1);
    end
    
    ht=triplot(TRI,S(1,:),S(2,:));
    harray(1)=ht(1);
    hf=get(h,'CurrentAxes');
    set(hf,'YDir','Reverse');
        set(hf,'FontSize',round(fontsize/2));
    set(ht,'Color',[1,0,0]);
    set(ht,'LineWidth',3);
    %title(methods{j},'Fontsize',fontsize);
if(j==1)
    hl=legend(harray,{'Ref. Shape','Trans. shapes'});
set(hl,'FontSize',fontsize,'Location','SouthWest');
end
hold off;
%print(h,'-depsc',sprintf('figs/exp2DFaceDatasetRef%s',methods{j}));
saveas(h,sprintf('figs/exp2DFaceDataset/exp2DFaceDatasetRef%s',methods{j}),'fig');
print(h,'-depsc',sprintf('figs/exp2DFaceDataset/exp2DFaceDatasetRef%s',methods{j}));
end

%% Generate tex report on results
filename='exp2DFaceDataset.tex'
fd=fopen(filename,'w');
cadena='\\begin{tabular}{cccc}\n'
fprintf(fd,cadena);
cadena='Method&RMS error& Time\\\\\\hline\n'
fprintf(fd,cadena);
for i=[1:nmethods]
    fprintf(fd,'%s&%f&%f \\\\\\hline\n',methods{i},sqrt(Result{i}.errorvec(end)),Result{i}.time(1))
end
cadena='\\end{tabular}\n'
fprintf(fd,cadena);
fclose(fd);
