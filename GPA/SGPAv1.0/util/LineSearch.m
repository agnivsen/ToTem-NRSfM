function [gammat]=LineSearch(X,D,DeltaX,V)
  %Coefficients of a 4th order polynomial
  %a0+a1*gamma+a2*gamma^2+a3*gamma^3=0
  n=length(D);
  [d,m]=size(D{1});
  [S,A,a]=unvectorize(X,d,m,n);
  [DeltaS,DeltaA,Deltaa]=unvectorize(DeltaX,d,m,n);
  a0=0;a2=0;
  a1=0;a3=0;
  for i=1:n
      for j=1:m
        if(V{i}(j,j))
        rij=D{i}(:,j)-A{i}*S(:,j)-a{i};
        focij=A{i}*DeltaS(:,j)+DeltaA{i}*S(:,j)+Deltaa{i};
        socij=DeltaA{i}*DeltaS(:,j);
        a0=a0-rij'*focij;
        a1=a1+norm(focij)^2-2*rij'*socij;
        a2=a2+3*focij'*socij;
        a3=a3+2*norm(socij)^2;
        end        
      end
  end
  
  gamma_roots=roots([a3,a2,a1,a0]);
  if(isreal(gamma_roots)) %only real roots--> Test the best !!
      X2=X+gamma_roots(1).*DeltaX;
      [S1,A1,a1]=unvectorize(X2,d,m,n);
      T1.S=S1;T1.A=A1;T1.a=a1;      
      error1=dataspaceerror(T1,D,V);
      X2=X+gamma_roots(2).*DeltaX;
      [S1,A1,a1]=unvectorize(X2,d,m,n);
      T1.S=S1;T1.A=A1;T1.a=a1;      
      error2=dataspaceerror(T1,D,V);
      X2=X+gamma_roots(3).*DeltaX;
      [S1,A1,a1]=unvectorize(X2,d,m,n);
      T1.S=S1;T1.A=A1;T1.a=a1;      
      error3=dataspaceerror(T1,D,V);
      [min_error,i]=min([error1,error2,error3]);
      gammat=gamma_roots(i);
  else %complex conjugate roots and a real root
    gammat=gamma_roots(find(imag(gamma_roots)==0));    
  end
  