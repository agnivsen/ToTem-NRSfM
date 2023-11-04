% AffineRegistration --- Generalized Procrustes Analysis with affine transformations..
%
% T=AffineRegistration(D,V)
%
% Arguments:
% 
% This function outputs the solution of the affine version of the GPA
% T is a structure with fields:
%       T.A is a cell-array with the rotational part of the transformations
%       T.a is a cell-array with the translational part of the
%       transformations
%       T.S is the reference shape.
%       T.errorvec is the data-space error obtained across iterations
%       T.niter is the number of iterations required.
%       T.time is the processing time of the algorithm.
%           T.time(1) is the total ellapsed time.
%           T.time(2) is the ellapsed time in the Initialization.
%           T.time(3) is the average time per iteration.
%           T.time(4) is the total time ellapsed in the refinement.
% D is the shape cell array D{i} is the dxm shape matrix of shape i.
% V is the matrix with missing data V(i,j)=1 if point j at shape i is missing. 
% options is an input structure with the following fields.
%      options.Ref={'ni','lm'};
%           'ni' => no iterative refinement.
%           'lm' => iterative refinement.
%      options.Init={'fct','ref'};
%           'fct' => Matrix factorization (no missing data)
%           'ref' => Initialization using reference space close form (DEFAULT).
%      options.method={'stratified','stratified'};
%       'alternation' => classical alternation method (only d=2,3))
%       'stratified'  => stratified approach. (DEFAULT)
%      options.verbose={'on','off'}
%       'on' => show debug information and error per iteration
%       'off'=> don't show debug information
%
%Copyright (C)2010  Adrien Bartoli (1), Daniel Pizarro (2) and Marco Loog
% (3)
% (1) Clermont Université - Clermont-Ferrand, France.
% (2) Universidad de Alcala- Alcala de Henares, Spain
% (3) Delft University - Nederlands
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

function T=AffineRegistration(D,V,options)

if(nargin==1)
    [D,V,options,n,d,m]=ProcessArgsAff(D);
elseif(nargin==2)
    [D,V,options,n,d,m]=ProcessArgsAff(D,V);
else
    [D,V,options,n,d,m]=ProcessArgsAff(D,V,options);
end
time=zeros(1,4);
tic;
switch(lower(options.Init))

    case 'ref' % Reference space method (default and allows missing data)
        KBi=cell(1,n);
        KBhat=cell(1,n);
        for i=1:n
            Di = [D{i} .* V{i}; ones(1,m)]; 
            Dih = Di';
            Dihat=Dih'*Dih;
            if(rcond(Dihat)<1e-12)
                Dihati=pinv(Dihat);
            else
                 Dihati=inv(Dihat);
            end
            KBi{i}=Dihati*Dih';
            KBhat{i}=Dih*KBi{i};
        end

        W=ones(m,m);
        for i=1:n
            k1=1;
            for j1=1:m
                k2=1;
                if(V{i}(j1)==1)
                    for j2=1:m
                        if(V{i}(j2)==1)
                            if(j1==j2)
                                W(j1,j2)=W(j1,j2)+1;
                            end
                            W(j1,j2)=W(j1,j2)-KBhat{i}(k1,k2);
                            k2=k2+1;
                        end

                    end

                    k1=k1+1;
                end
            end
        end
        
        W=W+norm(W, 'fro') * ones(m,m);
        % Extract the reference shape and the Affine transformations
        [U,Diag,Ut]=svd(W);
        S=U(:,size(U,2)-d+1:size(U,2));
        B=cell(1,n);
        for i=1:n
            B{i}=zeros(d+1,d);
            k=1;
            for j=1:m
                if(V{i}(j)==1)
                    for l1=1:d+1
                        for l2=1:d
                            B{i}(l1,l2)=B{i}(l1,l2)-KBi{i}(l1,k)*S(j,l2);
                        end
                    end
                    k=k+1;
                end
            end
            B{i}=B{i}';
        end

        %Extract the forward transformations
        A=cell(1,n);
        a=cell(1,n);
        for i=1:n
            A{i}=-inv(B{i}(1:d,1:d));
            a{i}=A{i}*B{i}(:,d+1);
        end
        S=S';

    case 'fct' %Affine factorization method (does not allow missing data)

        mu=zeros(d,n);
        for i=1:n
            for j=1:m
                mu(:,i)=mu(:,i)+D{i}(:,j);
            end
            mu(:,i)=mu(:,i)./m;
        end
        M=zeros(d*n,m);
        for j=1:m
            for i=1:n
                M(d*(i-1)+1:d*i,j)=D{i}(:,j)-mu(:,i);
            end
        end
        [Us,Ds,Vs]=svd(M);
        Ds=sqrt(Ds(1:d,1:d));
        Us=Us(:,1:d)*Ds;
        Vs=Vs';
        Vs=Ds*Vs(1:d,:);
        a=cell(1,n);
        A=cell(1,n);
        for i=1:n
            A{i}=Us(1+d*(i-1):i*d,1:d);
            a{i}=mu(:,i);
        end
        S=Vs;
end

% Iterative refinement methods in the data-space cost
% sum_i ||D_i-A_i S-a_i||^2_F

%Parameters X=[A{1}(:)',a{1}(:)',\cdots,A{n}(:)',a{n}(:)', S(:)];
%Initialize them using the reference-space cost solution
X=zeros(n*d*(d+1),1);
for i=1:n
    M=A{i}';
    X((d*d+d)*(i-1)+1:i*d*(d+1))=[M(:);a{i}];
end
X=[X;S(:)];
time(2)=toc;
switch lower(options.Ref)
    case 'ni'
        maxiter=1000;
        T.A=A;
        T.a=a;
        T.S=S;
        error=dataspaceerror(T,D,V);
        niter=1;
        errorvec=error.*ones(1,maxiter);
        switch(options.verbose)
            case 'on'
                switch(options.Init)
                    case 'fct'
                        disp(sprintf('Affine Factorization method  error=%f',error))
                    case 'ref'
                        disp(sprintf('Reference Space method  error=%f',error))
                end
        end
        
    case 'lm'
        %Levenberg-Mardquardt
        T.A=A;
        T.a=a;
        T.S=S;
        error=dataspaceerror(T,D,V);
        error_ant=0;
        maxiter=options.maxiter;
        niter=1;
        errorvec=zeros(1,maxiter);
        lambda=0.01;
        k=1;
        dispv(sprintf('LM  error=%f -- iteration=%d -- lambda=%f',error,niter,lambda),options.verbose);
        while(abs(error-error_ant)>10^(-6) && niter<maxiter)
            error_ant=error;
            errorvec(k)=error;
            [DeltaX,X2]=updateNormEqs(X,D,lambda,V);
            [S,A,a]=unvectorize(X2,d,m,n);T1.S=S;T1.A=A;T1.a=a;
            error2=dataspaceerror(T1,D,V);
            if(error2>error)
                lambda=lambda*10;
            else
                lambda=lambda/10;
            end
            niter=niter+1;
            X=X2;
            error=error2;
            dispv(sprintf('LM  error=%f -- iteration=%d -- lambda=%f',error,niter,lambda),options.verbose);
            k=k+1;
        end
        if(k>1)
            errorvec(k:end)=errorvec(k-1);
        end
end
%Set output solutions
[S,A,a]=unvectorize(X,d,m,n);
T.S=S;
T.A=A;
T.a=a;
T.errorvec=errorvec;
T.niter=niter;
time(1)=toc;
time(3)=toc-time(2);
time(4)=time(3)/niter;
T.time=time;
end



function [D,V,options,n,d,m]=ProcessArgsAff(D,V,options)

n=length(D);
[d,m]=size(D{1});
if(nargin<3)
    options.Ref='lm';
    options.Init='ref';
    options.verbose='off';
    options.maxiter=200;
    if(nargin<2)
        for i=1:n
            V{i}=ones(1,m);
        end
    end
else
    if(~iscell(V))
    V2=V;
    V=cell(1,n);
    for i=[1:n]
    V{i}=V2(i,:);
    end
    
    
    
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
    if(~isfield(options,'Init'))
        options.Init='ref';
    end

end
Vs=zeros(1, m);
for i=1:n
    Vs=Vs+V{i};
end
if(sum(Vs,2)==m*n)
 if(~strcmp(options.Ref,'ni'))
    dispv('No missing data...Switching to factorization method',options.verbose);
    options.Init='fct'; %modified by fang.bai
 end
else
    if(strcmp(options.Init,'fct'))
        disp('Error: Factorization method does not allow missing data.');
        disp('Error: Switching to reference-space method.');
        options.Init='ref';
    end
    
end


for i = 1:n
    [d1, d2] = size(V{i});
    if(d1 == d2)
        V{i} = diag(V{i})';
    end
end


end

