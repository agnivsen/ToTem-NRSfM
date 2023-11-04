function [length] = arcLength(p1, p2)
    % This computation is based on https://drive.google.com/file/d/1lRQ9y-0_f1_Co12m42DV0Q0KCnpYVk5f/view?usp=sharing
    % which is based on standard arc length computation on sphere
    % (e.g.: follow any diff. geom. text book, or this discussion: https://math.stackexchange.com/questions/231221/great-arc-distance-between-two-points-on-a-unit-sphere)

    assert(numel(p1)==2, 'Each input point must be a [2 x 1] array');
    assert(numel(p2)==2, 'Each input point must be a [2 x 1] array');

    theta1_1 = p1(1); theta1_2 = p2(1);
    theta2_1 = p1(2); theta2_2 = p2(2);

%     close all; pause(0.1);
%     sphere; hold on;
%     axis equal;
%     P1 = [-(cos(theta1_1).*cos(theta2_1)).' -sin(theta1_1).' -(cos(theta1_1).*sin(theta2_1)).'];
%     P2 = [-(cos(theta1_2).*cos(theta2_2)).' -sin(theta1_2).' -(cos(theta1_2).*sin(theta2_2)).'];
%     plot3(P1(1), P1(2), P1(3),'ro');
%     plot3(P2(1), P2(2), P2(3),'go');

    dotProd = (cos(theta1_1)*cos(theta1_2)*cos(theta2_1 - theta2_2)) + (sin(theta1_1)*sin(theta1_2));
    length = acos(dotProd);
end