% distance score between two shapes
% used to evaluate rigidity of transformations
function distScore = DistScore (D, Dtrans)

DT = delaunayTriangulation (D');

edges = DT.edges();

diff = zeros (size(edges,1), 1);

for jj = 1 : size(edges, 1)
    
    e = edges(jj, :);
    
    distD = norm( D(:, e(1)) - D(:, e(2)) );
    
    distDtrans = norm ( Dtrans(:, e(1)) - Dtrans(:, e(2)) );
    
    diff(jj) =  abs ((distD - distDtrans) / distD);
    
end

distScore = mean (diff);

end