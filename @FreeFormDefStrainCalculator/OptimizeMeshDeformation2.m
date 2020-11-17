%-------------------------------------------------------------------
% OptimizeDeformation:
% Iteratively fits the mesh to frames from MinFrame to MaxFrames.
% This is done by searching for the nodal parameters that best
% explains the displacement fields.
%
% OUTPUT:
% Xg, Yg : Fitted nodal parameters (num parameters x num frames)
%-------------------------------------------------------------------
function [Px, Py] = OptimizeMeshDeformation2 (this, NumFrames, Px0, Py0, DoOptimize, Iterations)
    if isempty(Px0) || isempty(Py0)
        % Initial guesses of the nodal parameters.
        % Number of nodes x Number of frames.
        Px(:,1) = this.Px_initial;
        Py(:,1) = this.Py_initial;
    else
        Px(:,1) = Px0;
        Py(:,1) = Py0;
    end
    
    [MeshX, MeshY] = meshgrid(0:(this.n * 5), 0:(this.n * 5));
    MaterialX = normalize(MeshX, 2, 'range') .* this.m;
    MaterialY = normalize(MeshY, 1, 'range') .* this.m;
    MaterialPoints = [MaterialX(:),  MaterialY(:)];

    [E, F] = this.GenerateConversionMatrices(MaterialPoints);
    NN = F';
    
    EulerPointsX = NN * Px(:,1);
    EulerPointsY = NN * Py(:,1);
    
    PointsArray = [EulerPointsX, EulerPointsY];
    
    for CurrentFrame = 2:NumFrames
        disp(['-- Processing frame ', num2str(CurrentFrame)]);
        
        Fx = this.DisplacementFields{CurrentFrame - 1,1}{1};
        Fy = this.DisplacementFields{CurrentFrame - 1,2}{1};
        
        PointsArray(:,1) = PointsArray(:,1) + Fx(PointsArray(:,1), PointsArray(:,2));
        PointsArray(:,2) = PointsArray(:,2) + Fy(PointsArray(:,1), PointsArray(:,2));
        
        fx = F * PointsArray(:,1);
        fy = F * PointsArray(:,2);
        
        Px(:,CurrentFrame) = inv(E' * E) * E' * fx;
        Py(:,CurrentFrame) = inv(E' * E) * E' * fy;
    end
end
