function hFig = VisualizeCameraAndPointclouds3D (fid, PoseCell, Points3D)

hFig = figure(fid);

if (nargin == 3)
    scatter3(Points3D(1,:), Points3D(2,:), Points3D(3,:), '.');
end

hold on;

TrajectoryArray = zeros(3, length(PoseCell));

for ii = 1 : length(PoseCell)
    
    % In our definition
    % R'*(Map - t) = Data
    % which means:    R*Data + t = Map
    % The pose is [R, t] here.
    P = PoseCell{ii};
    
    % For camera coordinate system
    % P = [R,  -R*C]
    % R*(Map - C) = Data
    % R' * Data + C = Map
    R = P(1:3, 1:3)';
    C = P(1:3, 4);
    
    TrajectoryArray(:,ii) = C;
    
    plotCamera('Location', C, 'Orientation', R, 'Opacity', 0, 'Color', 'b', 'AxesVisible', false, 'Size', 1.0);
    
end

plot3(TrajectoryArray(1,:), TrajectoryArray(2,:), TrajectoryArray(3,:), 'r-', 'LineWidth', 2);

hold off;
axis equal;
view(3);

pause(0.1);

end