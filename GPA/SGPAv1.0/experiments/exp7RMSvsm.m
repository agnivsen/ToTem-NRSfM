% Affine experiments:
addpath('../')
addpath('../util');
options.title='exp7RMSvsm';
options.transformations='EUC';
options.n=5;
options.d=3;
options.m=50;
options.tau=0;
options.sigma=sqrt(0.01);
options.varpar='m';
options.parvalues=round(linspace(4,50,10));
options.iterations=100;
options.methods={'AFF-FCT','AFF-REF','AFF-ALL'};
genExperiment(options);



