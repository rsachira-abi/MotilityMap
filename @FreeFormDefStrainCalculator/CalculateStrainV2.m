%-------------------------------------------------------------------
% CalculateStrain:
% Calculate the strain tensor for the mesh defined by the global
% mesh parameters Xg, Yg.
%
% INPUT:
% Xg, Yg: Global mesh parameters.
%
% OUTPUT:
% X, Y         : Global (euler) coordinates of the strain tensor.
% StrainTensor : Strain tensor for the global (euler) coordinate
%                location in X, Y.
%-------------------------------------------------------------------
function [X, Y, StrainTensor] = CalculateStrainV2 (this, Xg, Yg)
    X_element_initial = this.A * this.Xg_initial;
    Y_element_initial = this.A * this.Yg_initial;
    X_element = this.A * Xg;
    Y_element = this.A * Yg;
    
    [MaterialX, MaterialY] = meshgrid(linspace(0,1, 10), linspace(0,1, 10));
    MaterialX = MaterialX(:);
    MaterialY = MaterialY(:);
    
    k = 1;
    for e = 1:size(this.Elements,1)
        Ux_initial = X_element_initial(this.NumElementParameters*(e - 1) + 1:this.NumElementParameters*(e - 1) + this.NumElementParameters);
        Uy_initial = Y_element_initial(this.NumElementParameters*(e - 1) + 1:this.NumElementParameters*(e - 1) + this.NumElementParameters);
        Ux = X_element(this.NumElementParameters*(e - 1) + 1:this.NumElementParameters*(e - 1) + this.NumElementParameters);
        Uy = Y_element(this.NumElementParameters*(e - 1) + 1:this.NumElementParameters*(e - 1) + this.NumElementParameters);
        
        Points = this.CalculatePoint(MaterialX, MaterialY, Ux, Uy);
        PointsOld = this.CalculatePoint(MaterialX, MaterialY, Ux_initial, Uy_initial);
        
        dX_dx = zeros(size(MaterialX));
        dY_dx = zeros(size(MaterialX));
        dX_dy = zeros(size(MaterialX));
        dY_dy = zeros(size(MaterialX));
        
        dX_dx_undef = zeros(size(MaterialX));
        dY_dx_undef = zeros(size(MaterialX));
        dX_dy_undef = zeros(size(MaterialX));
        dY_dy_undef = zeros(size(MaterialX));
        
        for n = 1:this.NumElementParameters
            
            dX_dx = dX_dx + this.dsi_dx{n}(MaterialX, MaterialY) * Ux(n);
            dY_dx = dY_dx + this.dsi_dx{n}(MaterialX, MaterialY) * Uy(n);
            dX_dy = dX_dy + this.dsi_dy{n}(MaterialX, MaterialY) * Ux(n);
            dY_dy = dY_dy + this.dsi_dy{n}(MaterialX, MaterialY) * Uy(n);
            
            % Undeformed
            dX_dx_undef = dX_dx_undef + this.dsi_dx{n}(MaterialX, MaterialY) * Ux_initial(n);
            dY_dx_undef = dY_dx_undef + this.dsi_dx{n}(MaterialX, MaterialY) * Uy_initial(n);
            dX_dy_undef = dX_dy_undef + this.dsi_dy{n}(MaterialX, MaterialY) * Ux_initial(n);
            dY_dy_undef = dY_dy_undef + this.dsi_dy{n}(MaterialX, MaterialY) * Uy_initial(n);
        end
        
        E_11 = 0.5 .* ((dX_dx .* dX_dx) + (dY_dx .* dY_dx) - (dX_dx_undef .* dX_dx_undef) - (dY_dx_undef .* dY_dx_undef));
        E_12 = 0.5 .* ((dX_dx .* dX_dy) + (dY_dx .* dY_dy) - (dX_dx_undef .* dX_dy_undef) - (dY_dx_undef .* dY_dy_undef));
        E_21 = 0.5 .* ((dX_dy .* dX_dx) + (dY_dy .* dY_dx) - (dX_dy_undef .* dX_dx_undef) - (dY_dy_undef .* dY_dx_undef));
        E_22 = 0.5 .* ((dX_dy .* dX_dy) + (dY_dy .* dY_dy) - (dX_dy_undef .* dX_dy_undef) - (dY_dy_undef .* dY_dy_undef));

%         E_11 = 0.5 .* ((dX_dx .* dX_dx) + (dX_dy .* dX_dy) - (dX_dx_undef .* dX_dx_undef) - (dX_dy_undef .* dX_dy_undef));
%         E_12 = 0.5 .* ((dX_dx .* dY_dx) + (dX_dy .* dY_dy) - (dX_dx_undef .* dY_dx_undef) - (dX_dy_undef .* dY_dy_undef));
%         E_21 = 0.5 .* ((dY_dx .* dX_dx) + (dY_dy .* dX_dy) - (dY_dx_undef .* dX_dx_undef) - (dY_dy_undef .* dX_dy_undef));
%         E_22 = 0.5 .* ((dY_dx .* dY_dx) + (dY_dy .* dY_dy) - (dX_dx_undef .* dX_dx_undef) - (dY_dy_undef .* dY_dy_undef));
        
        a = (MaterialX == 1) & (MaterialY == 0); A = (MaterialX == 0) & (MaterialY == 0);
        b = (MaterialX == 0) & (MaterialY == 0); B = (MaterialX == 1) & (MaterialY == 1);
        c = (MaterialX == 0) & (MaterialY == 1); C = (MaterialX == 1) & (MaterialY == 0);
        d = (MaterialX == 0) & (MaterialY == 1); D = (MaterialX == 0) & (MaterialY == 0);
        
        dE_11 = hypot(PointsOld(a,1) - PointsOld(A,1), PointsOld(a,2) - PointsOld(A,2));
        dE_12 = hypot(PointsOld(b,1) - PointsOld(B,1), PointsOld(b,2) - PointsOld(B,2));
        dE_21 = hypot(PointsOld(c,1) - PointsOld(C,1), PointsOld(c,2) - PointsOld(C,2));
        dE_22 = hypot(PointsOld(d,1) - PointsOld(D,1), PointsOld(d,2) - PointsOld(D,2));
        
        % Assemble strain tensors
        for i = 1:length(MaterialX)
            
            X(k) = Points(i,1);
            Y(k) = Points(i,2);
            
            StrainTensor{k} = round([E_11(i) ./ (dE_11^2), E_12(i) ./ (dE_12^2); E_21(i) ./ (dE_21^2), E_22(i) ./ (dE_22^2)], 6);
            
            k = k + 1;
        end
    end
end