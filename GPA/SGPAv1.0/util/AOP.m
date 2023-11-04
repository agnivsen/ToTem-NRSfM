% AOP --- Solve the Absolute Orientation Problem.
%
% [R,t,s] = AOP(P1,P2,similarity,method)
%
% Arguments:
%
% P1 and P2 are dxm matrices composed of m d-dimensional vectors d=[2,3];
% The function search for the Rotation matrix R, traslation vector t and
% scale s so that P2=sRP1+t*ones(1,m)
%
% similarity={'on','off'}
%      'on' => Compute the Absolute Orientation Problem with scales (DEFAULT)
%      'off' => Compute the Absolute Orientation Problem with (s=1)
% method={'horn','SOS'}
%      'horn' => Uses the closed-form solution proposed by Horn: (DEFAULT)
%                "Closed-form solution of absolute orientation using unit
%                quaternions", BKP Horn - Journal of the Optical Society of
%                America A, 1987"
%      'SOS' => Uses the Sum Of Squares theory to compute the orientation
%      while ensuring positive scale (s>0).

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

function [R,T,s]=AOP(P1,P2,similarity,method)

if(nargin<4)
    method='horn';
    if(nargin<3)
        similarity='off';
    end
end
[d,n]=size(P1);

if(d>3)
    disp('Error: This function is designed for shapes dxm with d=[1,2]');
    R=[];
    T=[];
    s=[];
    return;
end
%%%First obtain the gravity center of each set and produce a normalized
%%%vectors

meanP1=mean(P1,2)*ones(1,n); 
P1n=P1-meanP1;
meanP2=mean(P2,2)*ones(1,n); 
P2n=P2-meanP2;
switch(lower(method))
    case 'horn'
        switch(d)
            case 2
                M=zeros(1,2);
                for i=1:n

                    M=M+[P1n(1,i)*P2n(1,i)+P1n(2,i)*P2n(2,i),-P1n(2,i)*P2n(1,i)+P1n(1,i)*P2n(2,i)];
                end
                q=M./norm(M); %%%Largest eigenvalue;
                R=[q(1),-q(2);q(2),q(1)];
                switch(similarity)
                    case 'on'
                        P1np=R*P1n;
                        %s=((P1np(:)'*P2n(:))/(P1n(:)'*P1n(:))); %signed
                        %scale.
                        s=sqrt(((P2n(:)'*P2n(:))/(P1n(:)'*P1n(:))));
                    case 'off'
                        s=1;
                end
                % Offset
                T=-s*R*P1+P2;
                T=mean(T,2);
            case 3
                %%%Obtain Matrix Sum of Products

                M=zeros(3,3);
                for i=1:n
                    M=M+P1n(:,i)*P2n(:,i)';
                end
                %%produce N Matrix for the solution
                N=zeros(4,4);
                N(1,:)=[M(1,1)+M(2,2)+M(3,3),M(2,3)-M(3,2),M(3,1)-M(1,3),M(1,2)-M(2,1)];
                N(2,:)=[M(2,3)-M(3,2),M(1,1)-M(2,2)-M(3,3),M(1,2)+M(2,1),M(3,1)+M(1,3)];
                N(3,:)=[M(3,1)-M(1,3),M(1,2)+M(2,1),-M(1,1)+M(2,2)-M(3,3),M(2,3)+M(3,2)];
                N(4,:)=[M(1,2)-M(2,1),M(3,1)+M(1,3),M(2,3)+M(3,2),-M(1,1)-M(2,2)+M(3,3)];
                [U,D]=eig(N);
                q=U(:,4);
                %q=V(:,1); %%%Largest eigenvalue;
                q=[q(1);-q(2:4)];
                R=quat2mat(q); %%%P2=RP1+T
                R=R';
                % Scale
                switch(similarity)
                    case 'on'
                        P1np=R*P1n;
                        %s=((P1np(:)'*P2n(:))/(P1n(:)'*P1n(:)));
                        s=sqrt(((P2n(:)'*P2n(:))/(P1n(:)'*P1n(:))));
                    case 'off'
                        s=1;
                end
                % Offset
                T=-s*R*P1+P2;
                T=[mean(T(1,:));mean(T(2,:));mean(T(3,:))];
        end

    case 'sos'
        switch(d)
            case 3
                % symbolical variables for quaternion
                syms q1 q2 q3 q4;
                r=[q1^2+q2^2-q3^2-q4^2; 2*q2*q3+2*q1*q4;  2*q2*q4-2*q1*q3; 2*q2*q3-2*q1*q4; q1^2+q3^2-q2^2-q4^2; 2*q3*q4+2*q1*q2; 2*q2*q4+2*q1*q3; 2*q3*q4-2*q1*q2; q1^2+q4^2-q3^2-q2^2];
                M=zeros(3*n,9);B=zeros(3*n,1);
                for i=1:n
                    M(3*(i-1)+1:3*i,:)=[[P1n(:,i)',zeros(1,6)];[zeros(1,3),P1n(:,i)',zeros(1,3)];[zeros(1,6),P1n(:,i)']];
                    B(3*(i-1)+1:3*i)=P2n(:,i);
                end
                Mtt=M'*M; Mt=-2*B'*M;
                % polynomial cost function
                cost=transpose(r)*Mtt*r+Mt*r;
                switch(similarity)
                    case 'off'
                        s=1;
                        equalities=q1^2+q2^2+q3^2+q4^2-1;
                        inequalities=q1;
                        degree=4;
                        % findbound in sostools
                        [gam,vars,opt]=findbound(cost,inequalities,equalities,degree,0);
                        R=quat2mat(opt)';
                    case 'on'
                        equalities=[];
                        inequalities=q1;
                        degree=4;
                        %find constrained bound in sostools
                        [gam,vars,opt]=findbound(cost,inequalities,equalities,degree,0);
                        s=norm(opt)^2;
                        R=quat2mat(opt./sqrt(s))';
                end
            case 2
                syms q1 q2;
                r=[q1;-q2;q2;q1];
                cost=0;
                for i=1:n
                    M=[[P1n(:,i)',zeros(1,2)];[zeros(1,2),P1n(:,i)']];
                    vec=M*r-P2n(:,i);
                    cost=cost+transpose(vec)*vec;
                end
                switch(similarity)
                    case 'off'
                        s=1;
                        equalities=q1^2+q2^2-1;
                        inequalities=q1;
                        degree=4;
                        [gam,vars,opt]=findbound(cost,inequalities,equalities,degree,0);
                        R=quat2mat(opt);
                    case 'on'
                        equalities=[];
                        inequalities=[];
                        degree=4;
                        [gam,vars,opt]=findbound(cost,inequalities,equalities,degree,0);
                        s=norm(opt);
                        R=[opt(1),-opt(2);opt(2),opt(1)]./s;
                end
        end
        T=-s*R*P1+P2;
        T=mean(T,2);
end
