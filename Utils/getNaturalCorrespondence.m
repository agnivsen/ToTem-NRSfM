function [flattenedPoints] = getNaturalCorrespondence(pointsOnToTem, topology)
    if(topology == 0)
        flattenedPoints = pointsOnToTem(:,1:2);
    elseif(topology==1)
        [r, theta] = flattenedCoordinatesCylinder(pointsOnToTem(:,1), pointsOnToTem(:,2), pointsOnToTem(:,3));
        flattenedPoints = [r (theta)]; 
    elseif(topology==2)
        [phi, theta] = flattenedCoordinatesSphere(pointsOnToTem(:,1), pointsOnToTem(:,2), pointsOnToTem(:,3));
        flattenedPoints = [phi theta]; 
    else
        error('This function has been adapetd only for the planar, cylindrical or spherical topologies');
    end
end


function [r, theta] = flattenedCoordinatesCylinder(X, Y, Z)
    r = Y;
    theta = atanDef(Z,X);
    for ii = 1:size(theta,1)
        if(theta(ii) > pi)
            theta(ii) = (theta(ii) - (2*pi));
        end
    end
end


function [phi, theta] = flattenedCoordinatesSphere(X, Y, Z)
    theta = atanDef(X,Z);
    for ii = 1:size(theta,1)
            theta(ii) = (theta(ii) - pi);
    end    
    phi = -asin(Y);
end

function [theta] = atanDef(X,Y)
    theta = [];
    for ii = 1:size(X,1)
        val = customAtan(Y(ii), X(ii));
        theta = [theta; val];
    end
end
function v=customAtan(y,x)
    if nargin==1 
        x=1;
    end
    v=nan;
    if x>0
        v=atan(y/x);
    end
    if y>=0 & x<0
        v=pi+atan(y/x);
    end
    if y<0 & x<0
        v=-pi+atan(y/x);
    end
    if y>0 & x==0
        v=pi/2;
    end
    if y<0 & x==0
        v=-pi/2;
    end
    if v<0
        v=v+2*pi;
    end
end