% Affine experiments:
addpath('../')
addpath('../util');
options.title='exp4RMSvsn';
options.n=5;
options.d=3;
options.m=50;
options.tau=0.5;
options.sigma=sqrt(0.01);
options.varpar='n';
options.parvalues=round(linspace(3,50,10));
options.iterations=100;
options.methods={'SIM-UPG','SIM-ALL','SIM-ALT'};
genExperiment(options);



