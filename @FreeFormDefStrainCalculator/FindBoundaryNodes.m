%-------------------------------------------------------------------
% FindBoundaryNodes:
% Find the boundary nodes to determine elements, based on the
% boundary points given.
%
% INPUT:
% NumNodes         : Number of nodes to use for fitting.
% PointsX, PointsY : User selected boundary points.
% ZeroX, ZeroY     : Hardcoded coordinates of zero length.
%
% OUTPUT:
% BoundaryNodes : Nodes of the boundary.
%-------------------------------------------------------------------
function [BoundaryNodes, BoundaryPoints, BoundaryPointsMaterial] = FindBoundaryNodes (~, NumNodes, Points, NodeMaterial)
    X = Points(:,1);
    Y = Points(:,2);

    Length = 0;
    for i = 2:length(X)
        Length(i) = Length(i - 1) + hypot(X(i - 1) - X(i), Y(i - 1) - Y(i));
    end
    
    [Length_, ia, ~] = unique(Length);
    BoundaryPointsMaterial = reshape(Length_ ./ Length(end), [], 1);

    LengthVec = NodeMaterial .* Length_(end);

    try
        Xs = reshape(spline(Length_, X(ia), LengthVec), [], 1);
        Ys = reshape(spline(X(ia), Y(ia), Xs), [], 1);
    catch
        Ys = reshape(spline(Length_, Y(ia), LengthVec), [], 1);
        Xs = reshape(spline(Y(ia), X(ia), Ys), [], 1);
    end

    BoundaryNodes = [Xs, Ys];
    BoundaryPoints = Points(ia,:);
end
