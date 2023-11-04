% GENEXPERIMENT --- Generate a testing experiment using artificially generated data
% for any method inside inside gpa function
% 
% genExperiment(options,name)
%
% INPUTS
%   options => structure defining the experiment
%       options.methods => cell array with methods to test. It can be any of the following
%                          or a group of them.
%                           {'AFF-FCT','AFF-REF','AFF-ALL','SIM-UPG','SIM-ALL',
%                           'SIM-ALT','EUC-UPG','EUC-ALL','EUC-ALT'}
%                          e.g. To test AFF-FCT and EUC-ALL
%                          options.methods={AFF-FCT,EUC-ALL};
%       options.transformations => Kind of transformation used to generate
%                                 artificial data. It can be any of the
%                                 following:
%                                   {'EUC','SIM','AFF} for
%                                   Euclidean,Similarity or Affine
%                                   respectively
%       options.n => number of shapes
%       options.d => dimension of each vector inside the shapes
%       options.m => number of vectors inside each shape
%       options.tau => percentage of outliers inside the shapes;
%       options.sigma => std of additive noise added to shape-data to simulate
%                       noise/rigidity. sigma = sqrt(0.01) is consider noise.
%                       sigma=sqrt(0.1) is consider deformations.
%       options.varpar => varing parameter. It can be any of the following
%            {'d','n','m','tau','sigmav','sigma'};
%       options.parvalues => rank of values in which the parameter is
%                           evaluated. E.g. options.varpar='m'; options.parvalues=[1:10];
%       options.iterations => number of iterations for each value of
%                             options.parvalues.
%       options.title => name of the file where the experiments are saved.
%  name => name of the file that contains an uncomplete experiment.
%  
% OUTPUTS
% 
% This function has no console outputs. It writes the results in the file
% mats/title-year-month-day-hour-minute-second.mat
% To visualize the results see function plotExperiment.

% Copyright (C)2010  Adrien Bartoli (1), Daniel Pizarro (2) and Marco Loog
% (3)
% (1) Clermont Université - Clermont-Ferrand, France.
% (2) Universidad de Alcala- Alcala de Henares, Spain
% (3) Delft University - Nederlands
% Send bug reports and feedback to: dani.pizarro@gmail.com
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA




function genExperiment(options,name)

if(nargin<2)
    %options=ProcessArgs(options);
    if(length(options)==0)
        disp('Exiting...');
        return
    end
    
    timevec=clock;
    n=options.n;
    d=options.d;
    m=options.m;
    tau=options.tau;
    sigma=options.sigma;
    N=options.iterations;
    nvalues=length(options.parvalues);
    nmethods=length(options.methods);
    resultpervalue=cell(N,nmethods);
    results=cell(1,nvalues);
    title=options.title;
    opti.verbose='off';
    opti.maxiter=250;
    savstring=sprintf('save mats/%s-%d-%d-%d-%d-%d-%d',options.title,timevec(1),timevec(2),timevec(3),timevec(4),timevec(5),round(timevec(6)));
    j=1;
    k=1;
    valueindex=1;
    nvalues0=valueindex;
    nmethods0=j;
    N0=k;
    index0=1;
else
    j=1;
    k=1;
    valueindex=1;
    index=1;
    eval(sprintf('load mats/%s',name));
    nvalues0=valueindex;
    nmethods0=j;
    N0=k;
    index0=index;
    savstring=sprintf('save mats/%s',name);
end

index=index0;
for valueindex=nvalues0:nvalues
    varingparam=options.varpar;
    eval(sprintf('%s=%f;',varingparam,options.parvalues(valueindex)));
    for k=N0:N
        S0=RandomShape(d,m);
        for i=[1:n]
            switch(options.transformations)
                case 'EUC'
                    [A0{i},a0{i}]=RandomEuclideanTx(d,1,'off');
                case 'SIM'
                    [A0{i},a0{i}]=RandomEuclideanTx(d,1,'on');
                case 'AFF'
                    [A0{i},a0{i}]=RandomAffineTx(d,1);
            end
            D{i}=A0{i}*S0+[a0{i}]*ones(1,m)+sigma*randn(d,m);
        end
        
        V=generate_missingdata(tau,n,m,d);
        for j=nmethods0:nmethods
            opti.method=options.methods{j};
            T=gpa(D,V,opti);
            T.error=dataspaceerror(T,D,V);
            resultpervalue{k,j}=T;
            eval(savstring);
            if(index>1)
            fprintf([repmat(8,1,length(cadena)+1) '']);
            end
            cadena=sprintf('Process %.2f/100.00',100*index/(N*nvalues*nmethods));
            disp(cadena);
            index=index+1;
        end
        nmethods0=1;
        j=j+1;    
    end
    N0=1;
    k=k+1;
  results{valueindex}=resultpervalue;
end
valueindex=valueindex+1;
eval(savstring);
end

function options=ProcessArgs(options)
    if(~isfield(options,'title'))
        options.title='Default';
    end
    if(~isfield(options,'transformations'))
        options.transformations='SIM';
    end
    if(~isfield(options,'n'))
        options.n=5;
    end
    if(~isfield(options,'m'))
        options.m=50;
    end
    if(~isfield(options,'d'))
        options.d=3;
    end
    if(~isfield(options,'tau'))
        options.tau=50;
    end
    if(~isfield(options,'sigma'))
        options.sigma=sqrt(0.01);
    end
        if(~isfield(options,'iterations'))
        options.iterations=100;
    end
    if(~isfield(options,'varpar'))
       disp('Error: varpar field missing');
       options=[];
       return;
    end
    if(~isfield(options,'parvalues'))
       disp('Error: parvalues filed missing');
       options=[];
       return;
    end
    if(~isfield(options,'methods'))
       disp('Error: parvalues filed missing');
       options=[];
       return;
    end
end



