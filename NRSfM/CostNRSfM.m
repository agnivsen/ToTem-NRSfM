classdef CostNRSfM<handle
    %% Method for computing DCF based cost for NRSfM 
    properties(Access=private)
        W = 0.5;
    end
    
    methods(Access = public)
        
        function [obj] = CostNRSfM(W)
            obj.W = W;
        end

        function [cost] = costSigned(obj, alpha, beta, theta, d)
            cost = [-(obj.W)*alpha ((1 - obj.W).*( alpha.^2 + beta.^2 - d.^2 - (2.*alpha.*beta.*cos(theta)) )).^2]; 
        end
    end
    
end
