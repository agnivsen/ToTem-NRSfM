function [largestDiagonal] = getLargestDiagonal(P)
    P = P - mean(P);
    minX = min(P(1,:)); maxX = max(P(1,:));
    minY = min(P(2,:)); maxY = max(P(2,:));
    minZ = min(P(3,:)); maxZ = max(P(3,:));
    largestDiagonal = sqrt( (maxX - minX)^2 + (maxY - minY)^2 + (maxZ - minZ)^2 );
end