% Affine experiments:
addpath('../')
addpath('../util');
options.title='exp1RMSvsd';
options.transformations='AFF';
options.n=5;
options.d=3;
options.m=50;
options.tau=0;
options.sigma=sqrt(0.01);
options.varpar='d';
options.parvalues=[1:10];
options.iterations=100;
options.methods={'AFF-FCT','AFF-REF','AFF-ALL'};
genExperiment(options);



