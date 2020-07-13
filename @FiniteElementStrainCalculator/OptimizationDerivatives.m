
function der_ssd = OptimizationDerivatives (this, Ux0, Uy0, Ux1, Uy1, MaterialPoints, Fx, Fy)
    der_ssd = zeros(this.NumElementParameters, 2, 2);
    
    Point0 = this.CalculatePoint(MaterialPoints(:,1), MaterialPoints(:,2), Ux0, Uy0);
    Point1 = this.CalculatePoint(MaterialPoints(:,1), MaterialPoints(:,2), Ux1, Uy1);
    
    E1 = MaterialPoints(:,1);
    E2 = MaterialPoints(:,2);
    
    Gx = feval(Fx, Point0);
    Gy = feval(Fy, Point0);
    
    [dGx_x, dGx_y] = differentiate(Fx, Point0);
    [dGy_x, dGy_y] = differentiate(Fy, Point0);
    
    point = this.CalculatePoint(MaterialPoints(:,1), MaterialPoints(:,2), Ux0, Uy0);
    
    for n = 1:this.NumElementParameters
        % dX
        der_ssd(n,1,1) = sum(2 .* (Point0(:,1) + Gx - Point1(:,1)) .* (this.si{n}(E1,E2) + dGx_x .* this.si{n}(E1,E2)) + 2 .* (Point0(:,2) + Gy - Point1(:,2)) .* (dGy_x .* this.si{n}(E1,E2)));
        der_ssd(n,2,1) = sum(2 .* (Point0(:,1) + Gx - Point1(:,1)) .* (-this.si{n}(E1,E2)));
    
        % dY
        der_ssd(n,1,2) = sum(2 .* (Point0(:,2) + Gy - Point1(:,2)) .* (this.si{n}(E1,E2) + dGy_y .* this.si{n}(E1,E2)) + 2 .* (Point0(:,1) + Gx - Point1(:,1)) .* (dGx_y .* this.si{n}(E1,E2)));
        der_ssd(n,2,2) = sum(2 .* (Point0(:,2) + Gy - Point1(:,2)) .* (-this.si{n}(E1,E2)));
    end
end