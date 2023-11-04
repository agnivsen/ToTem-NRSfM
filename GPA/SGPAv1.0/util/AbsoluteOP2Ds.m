function [R,T,q,s]=AbsoluteOP2Ds(P1,P2)
  %%%First obtain the gravity center of each set and produce a normalized vectors
    n=size(P1,2);
    meanP1=[mean(P1(1,:)).*ones(1,n);mean(P1(2,:)).*ones(1,n)];
    
    P1n=P1-meanP1;
    
    meanP2=[mean(P2(1,:)).*ones(1,n);mean(P2(2,:)).*ones(1,n)];
    
    P2n=P2-meanP2;
    %%%Obtain Matrix Sum of Products
    
    M=zeros(1,2);
    
    for i=[1:n]
    
        M=M+[P1n(1,i)*P2n(1,i)+P1n(2,i)*P2n(2,i),-P1n(2,i)*P2n(1,i)+P1n(1,i)*P2n(2,i)];
        %M=M+[P1n(1,i)*P2n(1,i),-P1n(2,i)*P2n(1,i);P1n(1,i)*P2n(2,i),P1n(2,i)*P2n(2,i)];
    end
      
    %%produce N Matrix for the solutio    
    %[U,D,V]=svd(M);
    q=M./norm(M); %%%Largest eigenvalue;
    R=[q(1),-q(2);q(2),q(1)];
      % Scale
    s=1;
    P1np=R*P1n;
%    s=((P1np(:)'*P2n(:))/(P1n(:)'*P1n(:)));
    s=sqrt(((P2n(:)'*P2n(:))/(P1n(:)'*P1n(:))));
    %q(1:3)=-q(1:3);
    T=-s.*R*mean(P1')'+mean(P2')';
     
    