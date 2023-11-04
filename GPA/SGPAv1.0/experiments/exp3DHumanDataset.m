% Experiment with the human eva dataset.
clear all;
addpath('../');
addpath('../util');
%load ../../datasets/3DMocap/Shapes3D
%D={D{1:7}};
%V={V{1:7}};
%n=size(D,2);    
%m=size(D{1},2);
%d=3;
step=1;
nframes=200;
m=20;
n=7;
load ../../datasets/3DMocap/scriptsAlvaro/out/S1/Walking_3.mat
 Markers=pose_seq;
x=Markers(1:step:nframes*step,1:m,1);
    y=Markers(1:step:nframes*step,1:m,2);
    z=Markers(1:step:nframes*step,1:m,3);
    i=1;
    for k=1:10:n*10
        
        D{i}=([x(k,:);y(k,:);z(k,:)]);
        V{i}=eye(m);
        i=i+1;
    end




methods={'AFF-ALL','SIM-ALL','SIM-ALT','EUC-ALL','EUC-ALT'}
nmethods=length(methods);
Result=cell(1,nmethods);
for i=1:nmethods
options.method=methods{i};
options.verbose='on';
Result{i}=gpa(D,V,options);
end
%% Show shapes

%paleta=hsv(nmethods);
paleta=0.8.*[1,0,0;0,1,0;0,0,1;1,1,0;0,1,1;0,0,0;0,1,0];
linestyle={'-','-.','--','-','-.',':','--',':'};
harray=zeros(1,nmethods+1);
fontsize=30;
for i=[1:n]
    h=figure(i);
    Di=D{i};
    ht=renderskeleton(Di,[1,0,0],2,'--');
    hold on;
    hf=get(h,'CurrentAxes');
    set(hf,'YDir','Reverse');
    daspect([1 1 1])
%    set(ht,'Color',[0,0,0]);
%    set(ht,'LineWidth',2,'LineStyle','--');
    harray(1)=ht(2);
    namearray{1}=sprintf('Shape %d',i);
    for j=1:nmethods
        S=Result{j}.S;
        A=Result{j}.A;
        a=Result{j}.a;
        Di=A{i}*S+a{i}*ones(1,m);
        ht=renderskeleton(Di,[1,0,0],2,linestyle{j});
        %ht=showbody(Di,V{i},[1,0,0],2);
        hold on
        harray(j+1)=ht(2);
        namearray{j+1}=methods{j};
        hf=get(h,'CurrentAxes');
        set(hf,'YDir','Reverse');
        set(ht,'Color',paleta(j,:));
        set(ht,'LineWidth',2);
    end
   hl=legend(harray,namearray);
   set(hl,'FontSize',fontsize,'Location','Best');
hold off
%print(h,'-depsc',sprintf('figs/exp3DHumanDatasetShape%d',i));
saveas(h,sprintf('figs/exp3DHumanDatasetShape%d',i),'fig');

end

%% Plot reference shapes for all methods.

for j=[1:nmethods]
h=figure(j+n);
hold on;
hl=zeros(1,2);
[error,T2]=referencespaceerror(Result{j},D,V); 
S=T2.S;
    A=T2.A;
    a=T2.a;
    for i=1:n
    Dip=inv(A{i})*D{i}-inv(A{i})*a{i}*ones(1,m);
        ht=renderskeleton(Dip,[1,0,0],2);
  
    hl(1)=ht(1);
    hc=get(h,'CurrentAxes');
    set(hc,'YDir','Reverse');
    daspect([1 1 1]);
    set(ht,'Color',[0,0,1]);
    set(ht,'LineWidth',1);
    end
    ht=renderskeleton(S,[1,0,0],2);
    hl(2)=ht(1);
    hf=get(h,'CurrentAxes');
    set(hf,'YDir','Reverse');
    set(ht,'Color',[1,0,0]);
    set(ht,'LineWidth',3);
    hl=legend(hl,{'Ref. Shape','Trans. shapes'});
    set(hl,'FontSize',fontsize,'Location','Best');

hold off;
if(j==1)
    input('Move the figure');
    CameraPosition=get(hc,'CameraPosition');
    CameraTarget=get(hc,'CameraTarget');
    CameraUpVector=get(hc,'CameraUpVector');
    CameraViewAngle=get(hc,'CameraViewAngle');
    Vaxis=axis;
    ha=hc;
else
    prop=set(ha);
    prop=rmfield(prop,'Children');
    prop=rmfield(prop,'Title');
    prop=rmfield(prop,'XLabel');
    prop=rmfield(prop,'YLabel');
    prop=rmfield(prop,'ZLabel');
    prop=rmfield(prop,'Parent');
    propnames=fieldnames(prop);
    for k=[1:size(propnames,1)]
        set(hc,propnames{k},get(ha,propnames{k}));
    end


end

%print(h,'-depsc',sprintf('figs/exp3DHumanDatasetRef%s',methods{j}));
saveas(h,sprintf('figs/exp3DHumanDatasetRef%s',methods{j}),'fig');

end

%% Generate tex report on results
filename='exp3DHumanDataset.tex';
fd=fopen(filename,'w');
cadena='\\begin{tabular}{cccc}\n';
fprintf(fd,cadena);
cadena='Method&RMS error& Time\\\\\\hline\n';
fprintf(fd,cadena);
for i=[1:nmethods]
    fprintf(fd,'%s&%f&%f \\\\\\hline\n',methods{i},sqrt(Result{i}.errorvec(end)),Result{i}.time(1));
end
cadena='\\end{tabular}\n';
fprintf(fd,cadena);
fclose(fd);

