%----------------------------------------------------------------------------
%----------------------------------------------------------------------------
function [Px, Py] = LeastSquaresMeshFit (this, KnownPoints, MaterialPoints)
    
    [E, F] = this.GenerateConversionMatrices(MaterialPoints);
    
    fx = F * KnownPoints(:,1);
    fy = F * KnownPoints(:,2);

    Px = inv(E' * E) * E' * fx;
    Py = inv(E' * E) * E' * fy;
end