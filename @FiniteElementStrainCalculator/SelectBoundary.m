%-------------------------------------------------------------------
% SelectBoundary:
% Let the user select the boundary of the mesh. Used for initial
% fitting.
%
% INPUT:
% NumNodes      : Number of nodes to use for fitting.
% ZeroX, ZeroY  : Hardcoded coordinates of zero length.
%
% OUTPUT:
% BoundaryNodes : Nodes of the boundary.
%-------------------------------------------------------------------
function [BoundaryNodes] = SelectBoundary (NumNodes, ZeroX, ZeroY, LastX, LastY)
    
    X = [];
    Y = [];
    while (true)
        [x, y, button] = ginput(1);
        if (button ~= 1)
            break;
        end
        
        X = [X; x];
        Y = [Y; y];
        
        hold on
        scatter(x, y, 'r', 'filled');
    end
    
    if (nargin >= 3)
        X = [ZeroX; X];
        Y = [ZeroY; Y];
    end
    
    if (nargin == 5)
        X = [X; LastX];
        Y = [Y; LastY];
    end
    
    Length = 0;
    for i = 2:length(X)
        Length(i) = Length(i - 1) + hypot(X(i - 1) - X(i), Y(i - 1) - Y(i));
    end
    LengthVec = linspace(0, Length(end), NumNodes);
    
    try
        Xs = spline(Length, X, LengthVec);
        Ys = spline(X, Y, Xs);
    catch
        Ys = spline(Length, Y, LengthVec);
        Xs = spline(Y, X, Ys);
    end
    
    BoundaryNodes = [Xs', Ys'];
end