% Affine experiments:
addpath('../')
addpath('../util');
options.title='exp3rigidRMSvstau'
options.transformations='AFF';
options.n=16;
options.d=3;
options.m=70;
options.tau=0;
options.sigma=sqrt(0.01);
options.varpar='tau';
options.parvalues=linspace(0,0.7,10);
options.iterations=100;
options.methods={'AFF-REF','AFF-ALL'};
genExperiment(options);



