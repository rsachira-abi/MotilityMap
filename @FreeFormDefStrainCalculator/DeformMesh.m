%-------------------------------------------------------------------
% DeformMesh:
% Deforms the mesh according to the global coordinates of the
% tracking points at timestep 't'.
%
% INPUT:
% t : Time step, to which the mesh should be deformed.
%
% OUTPUT:
% Xg, Yg : New global mesh parameters defining the fitted mesh.
%-------------------------------------------------------------------
function [Xg, Yg] = DeformMesh (this, t)
    
    for e = 1:size(this.Elements,1)
        
        PointsX = this.TrackingPoints(this.AssignedElement == e,1,t);
        PointsY = this.TrackingPoints(this.AssignedElement == e,2,t);
    
        xs = this.MaterialPoints(this.AssignedElement == e,1);
        ys = this.MaterialPoints(this.AssignedElement == e,2);
        
        for i = 1:this.NumElementParameters
            for j = 1:this.NumElementParameters
            
                E(i,j) = 0;
            
                for d = 1:length(xs)
                    E(i,j) = E(i,j) + this.si{i}(xs(d), ys(d)) * this.si{j}(xs(d), ys(d));
                end
            end
            
            fx(i) = 0;
            fy(i) = 0;
            for d = 1:length(xs)
                fx(i) = fx(i) + this.si{i}(xs(d), ys(d)) * PointsX(d);
                fy(i) = fy(i) + this.si{i}(xs(d), ys(d)) * PointsY(d);
            end
        end
        
        % Solve the element
        X_solve = pinv(E) * fx';
        Y_solve = pinv(E) * fy';
        
        ElementsX(this.NumElementParameters*(e - 1) + 1:this.NumElementParameters*(e - 1) + this.NumElementParameters) = X_solve;
        ElementsY(this.NumElementParameters*(e - 1) + 1:this.NumElementParameters*(e - 1) + this.NumElementParameters) = Y_solve;
    end
    
    % Solve the global mesh
    Xg = inv(this.A' * this.A) * this.A' * ElementsX';
    Yg = inv(this.A' * this.A) * this.A' * ElementsY';
end