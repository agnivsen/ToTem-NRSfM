addpath('../util');
addpath('../unstable');
addpath('../');
nameexp1RMSvsd='exp1RMSvsd-2010-11-23-8-46-46.mat';
nameexp1RMSvsm='exp1RMSvsm-2011-1-17-14-53-57.mat';
nameexp1RMSvsn='exp1RMSvsn-2010-11-23-8-49-44.mat';
nameexp1RMSvssigma='exp1RMSvssigma-2011-1-21-15-34-42.mat';

nameexp2RMSvsd='exp2RMSvsd-2010-11-23-8-52-30.mat';
nameexp2RMSvsm='exp2RMSvsm-2010-11-23-8-51-55.mat';
nameexp2RMSvsn='exp2RMSvsn-2010-11-23-8-51-25.mat';
nameexp2RMSvssigma='exp2RMSvssigma-2010-11-23-8-51-8.mat';

nameexp3nonrigidRMSvstau='exp3nonrigidRMSvstau-2010-11-23-10-51-33.mat';
nameexp3rigidRMSvstau='exp3rigidRMSvstau-2010-11-23-10-52-0.mat';

nameexp4RMSvsd='exp4RMSvsd-2011-1-20-16-24-56.mat';
nameexp4RMSvsm='exp4RMSvsm-2011-1-20-16-25-12.mat';
nameexp4RMSvsn='exp4RMSvsn-2011-1-20-16-25-55.mat';
nameexp4RMSvssigma='exp4RMSvssigma-2011-1-21-15-34-25.mat';
nameexp4RMSvstau='exp4RMSvstau-2011-1-20-16-26-33.mat';

nameexp5RMSvsd='exp5RMSvsd-2011-1-20-16-24-27.mat';
nameexp5RMSvsm='exp5RMSvsm-2011-1-20-16-24-13.mat';
nameexp5RMSvsn='exp5RMSvsn-2011-1-20-16-23-59.mat';
nameexp5RMSvssigma='exp5RMSvssigma-2011-1-20-16-23-23.mat';
nameexp5RMSvstau='exp5RMSvstau-2011-1-20-16-23-0.mat';

nameexp6RMSvsd='exp6RMSvsd-2011-1-20-16-20-57.mat';
nameexp6RMSvsm='exp6RMSvsm-2011-1-20-16-21-30.mat';
nameexp6RMSvsn='exp6RMSvsn-2011-1-20-16-21-49.mat';
nameexp6RMSvssigma='exp6RMSvssigma-2011-1-20-16-22-13.mat';
nameexp6RMSvstau='exp6RMSvstau-2011-1-20-16-22-32.mat';

%%
markervector={'square','^','diamond','o','p','>','x'};
paleta=[1,0,0;0,1,0;0,0,1;1,1,0;0,1,1;0,0,0];
%%{AFF-FCT,AFF-REF,AFF-ALL,SIM-ALL,SIM-UPG,SIM-ALT}
%Results on experiment 1. affine, full-data, rigid (sigma2 = .01). vary (noise range), n, m and d
fprintf('\\begin{tabular}{cccccc}\n');
%fprintf('Experiment & Method & min time (sec) & avg. time (sec) & max. time (sec) \\\\ \n');
fprintf('Experiment & Method & varing $d$ (sec) & varing $m$ (sec) & varing $n$ (sec) & varing $\\sigma$ & varing $\\tau$ \\\\ \n');

timevsparams=[];
[options,nvalues,timevsparam,errorvsparam]=plotExperiment2(nameexp1RMSvsd,1,'exp1RMSvsd',markervector,paleta,1);
timevsparams=[timevsparams,timevsparam];
N=options.iterations;
ext1=mean(timevsparam(:,N*(nvalues-1):N*nvalues)')';
[options,nvalues,timevsparam,errorvsparam]=plotExperiment2(nameexp1RMSvsm,2,'exp1RMSvsm',markervector,paleta,0);
timevsparams=[timevsparams,timevsparam];
ext2=mean(timevsparam(:,N*(nvalues-1):N*nvalues)')';
[options,nvalues,timevsparam,errorvsparam]=plotExperiment2(nameexp1RMSvsn,3,'exp1RMSvsn',markervector,paleta,0);
timevsparams=[timevsparams,timevsparam];
ext3=mean(timevsparam(:,N*(nvalues-1):N*nvalues)')';
[options,nvalues,timevsparam,errorvsparam]=plotExperiment2(nameexp1RMSvssigma,4,'exp1RMSvssigma',markervector,paleta,0);
timevsparams=[timevsparams,timevsparam];
ext4=mean(timevsparam(:,N*(nvalues-1):N*nvalues)')';
for j=[1:length(options.methods)]
%fprintf('1 & %s & %.2e & %.2e & %.2e \\\\\n',options.methods{j},min(timevsparams(j,:)),mean(timevsparams(j,:)),max(timevsparams(j,:)));
fprintf('1 & %s & %.2e & %.2e & %.2e& %.2e& \\\\\n',options.methods{j},ext1(j),ext2(j),ext3(j),ext4(j));
end
%% 
%Results on 2. affine, full-data, rigid (sigma2 = .01). vary (noise range), n, m and d
timevsparams=[];
[options,nvalues,timevsparam,errorvsparam]=plotExperiment2(nameexp2RMSvsd,5,'exp2RMSvsd',markervector,paleta,0);
timevsparams=[timevsparams,timevsparam];
ext1=mean(timevsparam(:,N*(nvalues-1):N*nvalues)')';
[options,nvalues,timevsparam,errorvsparam]=plotExperiment2(nameexp2RMSvsm,6,'exp2RMSvsm',markervector,paleta,0);
timevsparams=[timevsparams,timevsparam];
ext2=mean(timevsparam(:,N*(nvalues-1):N*nvalues)')';
[options,nvalues,timevsparam,errorvsparam]=plotExperiment2(nameexp2RMSvsn,7,'exp2RMSvsn',markervector,paleta,0);
timevsparams=[timevsparams,timevsparam];
ext3=mean(timevsparam(:,N*(nvalues-1):N*nvalues)')';
[options,nvalues,timevsparam,errorvsparam]=plotExperiment2(nameexp2RMSvssigma,8,'exp2RMSvssigma',markervector,paleta,0);
timevsparams=[timevsparams,timevsparam];
ext4=mean(timevsparam(:,N*(nvalues-1):N*nvalues)')';
for j=[1:length(options.methods)]
%fprintf('2 & %s & %.2e & %.2e & %.2e \\\\\n',options.methods{j},min(timevsparams(j,:)),mean(timevsparams(j,:)),max(timevsparams(j,:)));
fprintf('2 & %s & %.2e & %.2e & %.2e& %.2e& \\\\\n',options.methods{j},ext1(j),ext2(j),ext3(j),ext4(j));
end
%% 
markervector2={'^','diamond','o','p','>','x'};
paleta2=paleta(2:3,:);
%results on exp 3. affine, rigid and nonrigid. vary missing-data
timevsparams=[];
[options,nvalues,timevsparam,errorvsparam]=plotExperiment2(nameexp3nonrigidRMSvstau,9,'exp3nonrigidRMSvstau',markervector2,paleta2,0,2);
timevsparams=[timevsparams,timevsparam(:,1:nvalues*N)];
ext1=mean(timevsparam(:,N*(nvalues-1):N*nvalues)')';
[options,nvalues,timevsparam,errorvsparam]=plotExperiment2(nameexp3rigidRMSvstau,10,'exp3rigidRMSvstau',markervector2,paleta2,0,2);
timevsparams=[timevsparams,timevsparam(:,1:nvalues*N)];
ext2=mean(timevsparam(:,N*(nvalues-1):N*nvalues)')';
for j=[1:length(options.methods)]
%fprintf('3 & %s & %.2e & %.2e & %.2e \\\\\n',options.methods{j},min(timevsparams(j,:)),mean(timevsparams(j,:)),max(timevsparams(j,:)));
fprintf('3.1 & %s &&&&& %.2e \\\\\n',options.methods{j},ext1(j));
fprintf('3.2 & %s &&&&& %.2e \\\\\n',options.methods{j},ext2(j));
end
%%
%Results on 4
timevsparams=[];
[options,nvalues,timevsparam,errorvsparam]=plotExperiment2(nameexp4RMSvsd,11,'exp4RMSvsd',markervector,paleta,1);
timevsparams=[timevsparams,timevsparam];
ext1=mean(timevsparam(:,N*(nvalues-1):N*nvalues)')';
[options,nvalues,timevsparam,errorvsparam]=plotExperiment2(nameexp4RMSvsm,12,'exp4RMSvsm',markervector,paleta,0);
timevsparams=[timevsparams,timevsparam];
ext2=mean(timevsparam(:,N*(nvalues-1):N*nvalues)')';
[options,nvalues,timevsparam,errorvsparam]=plotExperiment2(nameexp4RMSvsn,13,'exp4RMSvsn',markervector,paleta,0);
timevsparams=[timevsparams,timevsparam];
ext3=mean(timevsparam(:,N*(nvalues-1):N*nvalues)')';
[options,nvalues,timevsparam,errorvsparam]=plotExperiment2(nameexp4RMSvssigma,14,'exp4RMSvssigma',markervector,paleta,0);
timevsparams=[timevsparams,timevsparam];
ext4=mean(timevsparam(:,N*(nvalues-1):N*nvalues)')';
[options,nvalues,timevsparam,errorvsparam]=plotExperiment2(nameexp4RMSvstau,15,'exp4RMSvstau',markervector,paleta,0,2);
timevsparams=[timevsparams,timevsparam(:,1:nvalues*N)];
ext5=mean(timevsparam(:,N*(nvalues-1):N*nvalues)')';
for j=[1:length(options.methods)]
%fprintf('4 & %s & %.2e & %.2e & %.2e \\\\\n',options.methods{j},min(timevsparams(j,:)),mean(timevsparams(j,:)),max(timevsparams(j,:)));
fprintf('4 & %s & %.2e & %.2e & %.2e& %.2e& %.2e \\\\\n',options.methods{j},ext1(j),ext2(j),ext3(j),ext4(j),ext5(j));
end

%%
%Results on 5
timevsparams=[];
[options,nvalues,timevsparam,errorvsparam]=plotExperiment2(nameexp5RMSvsd,16,'exp5RMSvsd',markervector,paleta,1);
timevsparams=[timevsparams,timevsparam];
ext1=mean(timevsparam(:,N*(nvalues-1):N*nvalues)')';
[options,nvalues,timevsparam,errorvsparam]=plotExperiment2(nameexp5RMSvsm,17,'exp5RMSvsm',markervector,paleta,0);
timevsparams=[timevsparams,timevsparam];
ext2=mean(timevsparam(:,N*(nvalues-1):N*nvalues)')';
[options,nvalues,timevsparam,errorvsparam]=plotExperiment2(nameexp5RMSvsn,18,'exp5RMSvsn',markervector,paleta,0);
timevsparams=[timevsparams,timevsparam];
ext3=mean(timevsparam(:,N*(nvalues-1):N*nvalues)')';
[options,nvalues,timevsparam,errorvsparam]=plotExperiment2(nameexp5RMSvssigma,19,'exp5RMSvssigma',markervector,paleta,0);
timevsparams=[timevsparams,timevsparam];
ext4=mean(timevsparam(:,N*(nvalues-1):N*nvalues)')';
[options,nvalues,timevsparam,errorvsparam]=plotExperiment(nameexp5RMSvstau,20,'exp5RMSvstau',markervector,paleta,0,2);
timevsparams=[timevsparams,timevsparam(:,1:nvalues*N)];
ext5=mean(timevsparam(:,N*(nvalues-1):N*nvalues)')';
for j=[1:length(options.methods)]
%fprintf('5 & %s & %.2e & %.2e & %.2e \\\\\n',options.methods{j},min(timevsparams(j,:)),mean(timevsparams(j,:)),max(timevsparams(j,:)));
fprintf('5 & %s & %.2e & %.2e & %.2e& %.2e& %.2e \\\\\n',options.methods{j},ext1(j),ext2(j),ext3(j),ext4(j),ext5(j));
end
%%
%Results on experiment 6
timevsparams=[];
[options,nvalues,timevsparam,errorvsparam]=plotExperiment2(nameexp6RMSvsd,21,'exp6RMSvsd',markervector,paleta,1);
timevsparams=[timevsparams,timevsparam];
ext1=mean(timevsparam(:,N*(nvalues-1):N*nvalues)')';
[options,nvalues,timevsparam,errorvsparam]=plotExperiment2(nameexp6RMSvsm,22,'exp6RMSvsm',markervector,paleta,0);
timevsparams=[timevsparams,timevsparam];
ext2=mean(timevsparam(:,N*(nvalues-1):N*nvalues)')';
[options,nvalues,timevsparam,errorvsparam]=plotExperiment2(nameexp6RMSvsn,23,'exp6RMSvsn',markervector,paleta,0);
timevsparams=[timevsparams,timevsparam];
ext3=mean(timevsparam(:,N*(nvalues-1):N*nvalues)')';
[options,nvalues,timevsparam,errorvsparam]=plotExperiment2(nameexp6RMSvssigma,24,'exp6RMSvssigma',markervector,paleta,0);
timevsparams=[timevsparams,timevsparam];
ext4=mean(timevsparam(:,N*(nvalues-1):N*nvalues)')';
[options,nvalues,timevsparam,errorvsparam]=plotExperiment2(nameexp6RMSvstau,25,'exp6RMSvstau',markervector,paleta,0,2);
timevsparams=[timevsparams,timevsparam(:,1:nvalues*N)];
ext5=mean(timevsparam(:,N*(nvalues-1):N*nvalues)')';
for j=[1:length(options.methods)]
%fprintf('6 & %s & %.2e & %.2e & %.2e \\\\\n',options.methods{j},min(timevsparams(j,:)),mean(timevsparams(j,:)),max(timevsparams(j,:)));
fprintf('6 & %s & %.2e & %.2e & %.2e& %.2e& %.2e \\\\\n',options.methods{j},ext1(j),ext2(j),ext3(j),ext4(j),ext5(j));
end

fprintf('\\end{tabular}\n')
