%-------------------------------------------------------------------
% CalculateStrain:
% Calculate the strain tensor for the mesh defined by the control points
% Px, Py.
%
% INPUT:
% Xg, Yg: Global mesh parameters.
%
% OUTPUT:
% X, Y         : Global (euler) coordinates of the strain tensor.
% StrainTensor : Strain tensor for the global (euler) coordinate
%                location in X, Y.
%-------------------------------------------------------------------
function [E, lambda] = CalculateStrain (this, Px, Py, MaterialPoints)
    %E = zeros(2, 2, size(MaterialPoints, 1)); % Strain tensors E(11), E(12), E(21), E(22)
    E = zeros(size(MaterialPoints, 1), 4);
    
    [dx_1, dx_2, dy_1, dy_2] = this.Derivatives(Px, Py, MaterialPoints);
    [dX_1, dX_2, dY_1, dY_2] = this.Derivatives(this.Px_initial, this.Py_initial, MaterialPoints);
    
    E(:,1) = 0.5 .* ((dx_1 .* dx_1) + (dy_1 .* dy_1) - (dX_1 .* dX_1) - (dY_1 .* dY_1));
    E(:,2) = 0.5 .* ((dx_1 .* dx_2) + (dy_1 .* dy_2) - (dX_1 .* dX_2) - (dY_1 .* dY_2));
    E(:,3) = 0.5 .* ((dx_2 .* dx_1) + (dy_2 .* dy_1) - (dX_2 .* dX_1) - (dY_2 .* dY_1));
    E(:,4) = 0.5 .* ((dx_2 .* dx_2) + (dy_2 .* dy_2) - (dX_2 .* dX_2) - (dY_2 .* dY_2));
    
    lambda = zeros(size(E, 1), 2); %[E11, E22]
    
    G_1 = [dX_1, dY_1];
    G_2 = [dX_2, dY_2];
    
    for i = 1:size(G_1,1)
        G1(i,:) = linsolve([G_1(i,1), G_1(i,2); G_2(i,1), G_2(i,2)], [1; 0])';
        G2(i,:) = linsolve([G_1(i,1), G_1(i,2); G_2(i,1), G_2(i,2)], [0; 1])';
    end
    
    G11 = G1(:,1).^2 + G1(:,2).^2;
    G22 = G2(:,1).^2 + G2(:,2).^2;
    
    dE_1_sqrd = dX_1.^2 + dY_1.^2;
    dE_2_sqrd = dX_2.^2 + dY_2.^2;
    dE_1_sqrd = 1 ./ G11;
    dE_2_sqrd = 1 ./ G22;
    
    lambda(:,1) = sqrt( (2 .* E(:,1) ./ dE_1_sqrd) + 1) - 1;
    lambda(:,2) = sqrt( (2 .* E(:,4) ./ dE_2_sqrd) + 1) - 1;
    
    tf = (imag(lambda(:,1)) ~= 0) & (imag(lambda(:,2)) ~= 0);
    lambda(tf,:) = NaN;
    lambda = real(lambda);
end



