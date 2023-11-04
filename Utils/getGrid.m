function [grid] = getGrid(pcl, xRange, yRange, zRange)
    xR = [min(pcl(:,1)) max(pcl(:,1))];
    yR = [min(pcl(:,2)) max(pcl(:,2))];
    zR = [min(pcl(:,3)) max(pcl(:,3))];

    xG = linspace(xR(1), xR(2), xRange);
    yG = linspace(yR(1), yR(2), yRange);
    zG = linspace(zR(1), zR(2), zRange);

    [X, Y, Z] = meshgrid(xG, yG, zG);

    grid = [X(:) Y(:) Z(:)];
end