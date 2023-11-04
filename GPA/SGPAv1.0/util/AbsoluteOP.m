function [R,T,q]=AbsoluteOP(P1,P2)
  %%%First obtain the gravity center of each set and produce a normalized vectors
    n=size(P1,2);
    meanP1=[mean(P1(1,:)).*ones(1,n);mean(P1(2,:)).*ones(1,n);mean(P1(3,:)).*ones(1,n)];
    
    P1n=P1-meanP1;
    
        meanP2=[mean(P2(1,:)).*ones(1,n);mean(P2(2,:)).*ones(1,n);mean(P2(3,:)).*ones(1,n)];
    
    P2n=P2-meanP2;
    %%%Obtain Matrix Sum of Products
    
    M=zeros(3,3);
    
    for i=[1:n]
      M=M+P1n(:,i)*P2n(:,i)';
    end
      
    %%produce N Matrix for the solution
    N=zeros(4,4);
    N(1,:)=[M(1,1)+M(2,2)+M(3,3),M(2,3)-M(3,2),M(3,1)-M(1,3),M(1,2)-M(2,1)];
    N(2,:)=[M(2,3)-M(3,2),M(1,1)-M(2,2)-M(3,3),M(1,2)+M(2,1),M(3,1)+M(1,3)];
    N(3,:)=[M(3,1)-M(1,3),M(1,2)+M(2,1),-M(1,1)+M(2,2)-M(3,3),M(2,3)+M(3,2)];
    N(4,:)=[M(1,2)-M(2,1),M(3,1)+M(1,3),M(2,3)+M(3,2),-M(1,1)-M(2,2)+M(3,3)];
    
    [U,D,V]=svd(N);
    q=V(:,1); %%%Largest eigenvalue;
    q=[q(1);-q(2:4)];
    R=quat2mat(q);%real(quat2Rot(q)); %%%P2=RP1+T
    R=R';
    % Scale
    s=1;
    P1np=R*P1n;
    %s=P2n(:)'*(P1np(:))/(P1n(:)'*P1n(:));
    % Offset
    
    %q(1:3)=-q(1:3);
    T=-s*R*P1+P2;
    %[U1,D1,V1]=svd(T);
    %D1(2:end)=0;
    %T=U1*D1*V1';
    %T=U1*D1(:,1);
    T=[mean(T(1,:));mean(T(2,:));mean(T(3,:))];
    
    