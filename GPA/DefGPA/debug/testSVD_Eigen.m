clear all
close all
clc

DIM1 = 8;
DIM2 = 6;


B = 100 * rand(DIM1, DIM2); 

A = magic (DIM1);
A = (A+A')/2

[U, S, V] = svd (A)

fprintf(2, '[U, S, V] = svd (A), If A is sysmetric, then U^T * V = V^T * U. Note U ~= V');

U'*V - V'*U

[Q, E]= eig(A);
[ev, index] = sort (diag(E), 'descend');
Q = Q(:, index)

norm(A - Q * diag(ev) * Q')
norm(Q * Q' )


K = B' * B;
[Q, D] = eig ( (K+K')/2 );
[Ds, Is]= sort( diag (D), 'descend' );
Q = Q(:, Is);
D = diag( Ds )


[U, S, V] = svd (B);

D
SS = S .* S

K - Q * D * Q'

K - V * D * V'

Q * Q'
V * V'

abs (Q) - abs (V)
(Q + V) .* (Q-V)





[V1, E1, EA1] = byEigen (B' * B, 3);

[V2, E2, EA2] = bySVD (B, 3);

[E1, E2]

[sqrt(EA1), EA2]






clear all
close all
clc


R1 = SO3.Exp(10*randn(3,1));
D1 = diag (rand(3,1));

A = R1 * D1
B = D1 * R1






clear all
close all
clc

B = 10 * rand(3, 5); 

K = B * B';
[Q, D] = eig ( (K+K')/2 );

for i = 1 : size(B,2)
    pp{i} = B(:,i) * B(:,i)';
    cc{i} = Q' * pp{i} * Q;
end


X = zeros (3, 3);
for i = 1 : size(B, 2)
    X = X + cc{i}
end











        function [V, eigValEigen, eigenValueAll] = byEigen (K, d)

            [Q, D ] = eig( (K+K')/2 );
            
            eigVal = diag(D);
            
            [eigVal, Index] = sort(eigVal, 'descend');  Q = Q (:, Index);

            V = Q (:, end-d+1:end);
            
            eigenValueAll = eigVal;
            
            eigValEigen =  eigVal (end-d+1:end);

        end

        
        function [V, eigValSVD, singularValueAll] = bySVD (K, d)

            [~, D, V] = svd(K);
            
            V = V(:, end-d+1:end);
            
            singularValueAll = diag(D);
            
            eigValSVD = singularValueAll (end-d+1:end).^2;

        end
