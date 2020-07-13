%-------------------------------------------------------------------
% ElementBoundary:
% Output the set of points representing the boundary of the element
% described by Ux, and Uy.
%
% INPUT:
% Ux, Uy   : Element parameters that describe the element.
%
% OUTPUT:
% (Xs, Ys) : Set of coordinates (X, Y) of the boundary of the
%            element.
%-------------------------------------------------------------------
function [Xs, Ys]  = ElementBoundary (this, Ux, Uy)
    
    i = 1;
    for v = linspace(0, 1, 100)
        Point{1}(i,:) = this.CalculatePoint(v, 0, Ux, Uy);
        Point{2}(i,:) = this.CalculatePoint(0, v, Ux, Uy);
        Point{3}(i,:) = this.CalculatePoint(v, 1, Ux, Uy);
        Point{4}(i,:) = this.CalculatePoint(1, v, Ux, Uy);
        i = i + 1;
    end
    
    Xs = [Point{1}(:,1); Point{2}(:,1); flip(Point{3}(:,1)); flip(Point{4}(:,1))];
    Ys = [Point{1}(:,2); Point{2}(:,2); flip(Point{3}(:,2)); flip(Point{4}(:,2))];
end