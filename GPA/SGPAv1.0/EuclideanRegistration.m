% EUCLIDEANREGISTRATION --- Generalized Procrustes Analysis with similarity/euclidean transformations.
% This function outputs the solution of the rigid version of the GPA
%
% T=EuclideanRegistration(D,V)
%
% T is a structure with fields:
%       T.A is a cell-array with the rotational part of the transformations
%       T.a is a cell-array with the translational part
%       T.S is the reference shape.
%       T.errorvec is the data-space error obtained across iterations
%       T.niter is the number of iterations required.
%       T.time is the processing time of the algorithm.
%           T.time(1) is the total ellapsed time.
%           T.time(2) is the ellapsed time in the Initialization+upgrade from affine.
%           T.time(3) is the average time per iteration.
%           T.time(4) is the total time ellapsed in the refinement.
% D is the shape cell array D{i} is the dxm shape matrix of shape i with d=[1,2].
% V is the matrix with missing data V(i,j)=1 if point j at shape i is
% missing.
% options is an input structure with the following fields.
%      options.similarity={'on','off'};
%           'on' => Scales included in the problem
%           'off' => Scales not included. (DEFAULT)
%      options.Ref={'ni','lm'};
%           'ni' => no iterative refinement after update
%           'lm' => iterative refinement after update
%      options.upgrade={'fromreg','fromshape'};
%           'fromreg' => update from affine using registration (DEFAULT)
%           'fromshape' => update from affine using shapes.
%      options.method={'stratified','stratified'};
%       'alternation' => classical alternation method
%       'stratified'  => stratified approach. (DEFAULT)
%       options.maxiter => maximum number of iterations (=200 DEFAULT)

% Copyright (C)2010  Daniel Pizarro (1) and Adrien Bartoli (2)
% (1) Universidad de Alcala- Alcala de Henares, Spain
% (2) Clermont Université - Clermont-Ferrand, France.
%
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

function T=EuclideanRegistration(D,V,options)

% Check arguments
if(nargin==1)
    [D,V,options,n,d,m]=ProcessArgsEuc(D);
elseif(nargin==2)
    [D,V,options,n,d,m]=ProcessArgsEuc(D,V);
else
    [D,V,options,n,d,m]=ProcessArgsEuc(D,V,options);
end
if(d>3)
    disp('Error: This function is designed for shapes dxm with d=[1,2]');
    T.A=[];
    T.a=[];
    T.S=[];
    T.errorvec=[];
    T.niter=0;
    T.time=[];
    return;
end
time=zeros(1,4);
tic;
% Some useful constants
beta=zeros(1,n);
for i=1:n
    beta(i)=sum(V{i});
end

switch(options.method)

    case 'stratified'
        opti.verbose=options.verbose;
        T=AffineRegistration(D,V,options.affine);
        S=T.S;
        A=T.A;
        a=T.a;
        %% Upgrading to Euclidean
        switch(options.upgrade)
            case 'fromshape' %Upgrading from shape.
                %first we build the centred measurement matrix XM.
                XM=zeros(d*n,m);
                for j=1:m
                    for i=1:n
                        if(V{i}(j))
                            XM(d*(i-1)+1:d*i,j)=D{i}(:,j)-T.a{i};
                        else
                            XM(d*(i-1)+1:d*i,j)=T.A{i}*S(:,j); % Predicted measurements
                        end
                    end
                end
                Y=S*XM'*XM*S';
                Z=chol(Y);
                Q=XM*S'*inv(Z);
                for i=1:n
                    Qt=(sqrt(n)).*Q(d*(i-1)+1:i*d,:);
                    %   Qt=T.A{i}*inv(Z);
                    [Ut,Dt,Vt]=svd(Qt);
                    switch(options.similarity)
                        case 'off'
                            A{i}=Ut*Vt';
                        case 'on'
                            A{i}=Ut*Vt'*mean(diag(Dt));
                    end
                end
                S=Z*S./sqrt(n);
                T.S=S;
                T.A=A;
                T.a=a;

            case 'fromreg' %Upgrading from registration parameters (classical approach).

                %Test for consistent orientation
                A=T.A;
                a=T.a;
                S=T.S;
                w=sign(det(A{1}));
                % We divide the shapes into two groups. Those that agree
                % with w and those that don't agree.
                gindexes=ones(1,n);
                for i=2:n
                    if(w~=sign(det(A{i})))
                        gindexes(i)=0;
                        dispv(sprintf('[shape %d/%d]:not consistent signs',i,n),options.verbose);

                    end
                end
                %find which group is bigger.
                if(sum(gindexes)>sum(~gindexes))
                    gindexes=~gindexes;
                end
                %select the sign of the bigger group
                wp=find(~gindexes);
                w=sign(det(A{wp(1)}));
                %wp=[1:length(gindexes)];
                 
                switch(options.similarity)
                    case 'off'
                        M=(A{wp(1)})'*(A{wp(1)});
                        for i=2:length(wp)
                            M=M+(A{wp(i)})'*(A{wp(i)});
                        end
                        M=M./length(wp);
                    case 'on'
                        M=(A{wp(1)})'*(A{wp(1)});
                        M=M.*(1/det(M))^(1/d);
                        for i=2:length(wp)
                            Mtt=(A{wp(i)})'*(A{wp(i)});
                            M=M+Mtt.*(1/det(Mtt))^(1/d);
                        end
                end

                Z=chol(M);
                Z(d,d)=sign(det(Z))*w;
                Zi=inv(Z);

                %Applying the upgrade transformation
                E=cell(1,n);
                for i=1:length(wp)
                    Etemp=A{wp(i)}*(Zi)';
                    [Ue,De,Ve]=svd(Etemp);
                    switch(options.similarity)
                        case 'on'
                            E{wp(i)}=(1/d)*(trace(De)).*Ue*Ve';
                        case 'off'
                            E{wp(i)}=Ue*Ve';
                    end
                end
                T.a=a;
                T.A=E;
                A=E;
                S=Z'*S;
                T.S=S;

                % Recompute the euclidean registration parameters for the loosing
                % group
                for i=find(gindexes)
                    
                    %Compute visible part of the shape
                    meanshape =S;
                    [Rt,Tt,s]=AOP(meanshape(:, logical(V{i})), D{i}(:, logical(V{i})), options.similarity,'horn');
                    A{i}=s.*Rt;
                    a{i}=Tt;
                end
                %check again sign agreement.
                w=sign(det(A{1}));
                for i=2:n
                    if(w~=sign(det(A{i})))
                        disp(sprintf('[2ERROR shape %d/%d]:not consistent signs',i,n));
                    end
                end

                T.A=A;
                T.S=S;
                T.a=a;
        end
        time(2)=toc;
        switch lower(options.Ref)
            case 'ni' % no refinement
                maxiter=options.maxiter;
                error=dataspaceerror(T,D,V);
                dispv(sprintf('EUC-UPG  error=%f',error),options.verbose);
                niter=1;
                errorvec=error.*ones(1,maxiter);
                T.niter=niter;
                T.errorvec=errorvec;

            case 'lm' % Levenberg-Mardquardt iterations
                maxiter=options.maxiter;
                X=vectorizeEuc(S,A,a,d,m,n,options.similarity);
                T.A=A;T.a=a;T.S=S;
                error=dataspaceerror(T,D,V);
                error_ant=0;
                niter=1;
                errorvec=zeros(1,maxiter);
                lambda=0.01;
                k=1;
                dispv(sprintf('EUC-ALL  error=%f -- iteration=%d -- lambda=%f',error,niter,lambda),options.verbose);
                while(abs(error-error_ant)>10^(-6) && niter<maxiter)
                    error_ant=error;
                    errorvec(k)=error;
                    [DeltaX,X2]=updateNormEqsEuc(X,D,lambda,V,options.similarity);
                    [S,A,a]=unvectorizeEuc(X2,d,m,n,options.similarity);T1.S=S;T1.A=A;T1.a=a;
                    error=dataspaceerror(T1,D,V);
                    if(error>error_ant)
                        lambda=lambda*2;
                       
                    else
                        lambda=lambda/2;
                        X=X2;
                    end
                    niter=niter+1;
                    dispv(sprintf('EUC-ALL  error=%f -- iteration=%d -- lambda=%f',error,niter,lambda),options.verbose);
                    k=k+1;
                end
                if(k>1)
                    errorvec(k:end)=errorvec(k-1);
                end
                T.A=A;
                T.a=a;
                T.S=S;
                T.errorvec=errorvec;
                T.niter=niter;
        end

    case 'alternation' %Alternation approach


        incS=inf;
        itermax=options.maxiter;
        errorvec=zeros(1,itermax);
        meanshapeb=zeros(d,m);
        % Initilize all registration parameter to identity transformations
        A=cell(1,n);
        a=cell(1,n);
        for i=1:n
            A{i}=eye(d);
            a{i}=zeros(d,1);
        end
        iter=1;
        while(abs(incS)>10^(-6) && iter<itermax)

            %compute mean shape.
            meanshape=zeros(d,m);
            meanV=zeros(1,m);
            for i=1:n
                meanV=meanV+V{i};
                meanshape=meanshape+(inv(A{i})*D{i}-inv(A{i})*a{i}*ones(1,m)) .* V{i};
            end
            if(iter>0)
                meanshape=meanshape ./ meanV;
            else
                meanshape=inv(A{1})*D{1}-inv(A{1})*a{1}*ones(1,m);
            end
            % Increment of the reference shape
            incS=norm(meanshapeb-meanshape);
            %compute registration one-to-one in turns
            meanshapeb=meanshape;
            for i=1:n
                % Truncate the shape with only visible points
                [Rt,Tt,s]=AOP(meanshape(:, logical(V{i})), D{i}(:, logical(V{i})),options.similarity,'horn');
                A{i}=s.*Rt;
                a{i}=Tt;
            end

            T.A=A;
            T.a=a;
            T.S=meanshape;
            error=dataspaceerror(T,D,V);
            dispv(sprintf('ALT: error=%f, iter=%d, incS=%f',error,iter,incS),options.verbose);
            errorvec(iter)=error;
            iter=iter+1;
        end
        T.niter=iter;
        errorvec(iter:end)=errorvec(iter-1);
        T.errorvec=errorvec;

    case 'alternationref'

        incS=inf;
        itermax=options.maxiter;
        errorvec=zeros(1,itermax);
        meanshapeb=zeros(d,m);
        A=cell(1,n);
        a=cell(1,n);
        for i=1:n
            A{i}=eye(d);
            a{i}=zeros(d,1);
        end
        iter=1;
        while(abs(incS)>10^(-6) && iter<itermax)
            %compute mean shape.
            meanshape=zeros(d,m);
            meanV=zeros(1,m);
            for i=1:n
                meanV=meanV+V{i};
                meanshape=meanshape+(A{i}*D{i}+a{i}*ones(1,m)) .* V{i};
            end
            if(iter>1)
                meanshape=meanshape ./ meanV;
            else
                meanshape=A{1}*D{1}+a{1}*ones(1,m);
            end
            incS=norm(meanshapeb-meanshape);
            %compute registrataion one-to-one
            meanshapeb=meanshape;
            for i=1:n
                [Rt,Tt,s]=AOP(D{i}(:, logical(V{i})), meanshape(:, logical(V{i})),options.similarity,'horn');
                A{i}=s.*Rt;
                a{i}=Tt;
            end

            for i=1:n
                T.A{i}=inv(A{i});
                T.a{i}=-inv(A{i})*a{i};
            end
            T.S=meanshape;
            error=dataspaceerror(T,D,V);
            errorvec(iter)=error;
            dispv(sprintf('ALT: error=%f, iter=%d, incS=%f',error,iter,incS),options.verbose);
            iter=iter+1;
        end
        for i=1:n
            A{i}=inv(A{i});
            a{i}=-A{i}*a{i};
        end
        T.A=A;
        T.a=a;
        T.S=meanshape;
        errorvec(iter:end)=errorvec(iter-1);
        T.errorvec=errorvec;
        T.niter=iter;

end
time(1)=toc;
time(3)=toc-time(2);
time(4)=time(3)/T.niter;
T.time=time;
end

function [D,V,options,n,d,m]=ProcessArgsEuc(D,V,options)

n=length(D);
[d,m]=size(D{1});
if(nargin<3)
    options.method='stratified';
    options.Ref='lm';
    options.verbose='off';
    options.upgrade='fromreg';
    options.similarity='off';
    options.affine=[];
    if(nargin<2) %No missing data
        for i=1:n
            V{i}=ones(1,m);
        end
    end
else
    
    
        for i=1:n
            V{i}=ones(1,m);
        end
    
    if(~iscell(V))
        V2=V;
        V=cell(1,n);
        for i=[1:n]
            V{i}=diag(V2(i,:));
        end
    end
    
    if(~isfield(options,'method'))
        options.method='stratified';
    end
    if(~isfield(options,'Ref'))
        options.Ref='lm';
    end
    if(~isfield(options,'maxiter'))
        options.maxiter=200;
    end
    if(~isfield(options,'verbose'))
        options.verbose='off';
    end
    if(~isfield(options,'upgrade'))
        options.upgrade='fromreg';
    end
    if(~isfield(options,'similarity'))
        options.similarity='off';
    end
    if(~isfield(options,'affine'))
        options.affine.verbose=options.verbose;
    end
    
end

end





