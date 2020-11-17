
%-------------------------------------------------------------------
% OptimizeDeformation:
% Iteratively fits the mesh to frames from MinFrame to MaxFrames.
% This is done by searching for the nodal parameters that best
% explains the displacement fields.
%
% OUTPUT:
% Xg, Yg : Fitted nodal parameters (num parameters x num frames)
%-------------------------------------------------------------------
function [Px, Py] = OptimizeMeshDeformationCUDA (this, NumFrames, Px0, Py0, DoOptimize, Iterations)
    if isempty(Px0) || isempty(Py0)
        % Initial guesses of the nodal parameters.
        % Number of nodes x Number of frames.
        Px(:,1) = this.Px_initial;
        Py(:,1) = this.Py_initial;
    else
        Px(:,1) = Px0;
        Py(:,1) = Py0;
    end

    PieceSize = this.PIECE_SIZE;

    StepSize = floor(PieceSize / 2);
end