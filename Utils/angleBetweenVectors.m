function [angle] = angleBetweenVectors(u,v)
%% Find angle between two vectors u and v in radians
    angle = atan2(norm(cross(u,v)),dot(u,v));
end