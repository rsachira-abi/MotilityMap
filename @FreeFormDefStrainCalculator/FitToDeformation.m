%-------------------------------------------------------------------
% FitToDeformation:
% Iteratively fits the mesh to frames from MinFrame to MaxFrames.
% This is done by searching for the nodal parameters that best
% explains the displacement fields.
%
% INPUT:
% XgYg0   : Current guess of guess of X and Ynodal parameters
%           that describe the deformation.
%           (num parameters x num frames)
% MinFrame: Starting frame
% MaxFrame: Last frame
%
% OUTPUT:
% metric : Fitness of the current parameters.
%-------------------------------------------------------------------
function metric = FitToDeformation (this, XgYg0, MinFrame, MaxFrame)
        metric = 0;
        XgYg0 = reshape(XgYg, size(this.A, 2), [], 2);

        v = linspace(0, 1, 4);
        [MaterialX, MaterialY] = meshgrid(v, v);

        for i = 1:MaxFrame
                ElementX_i = this.A * XgYg(:,i,1);
                ElementY_i = this.A * XgYg(:,i,2);
                
                for j = 1:length(DisplacementFields{i,3})
                        if (DisplacementFields{i,4}{j} < this.ERROR_CUTOFF && DisplacementFields{i,3}{j} >= MinFrame && DisplacementFields{i,3}{j} <= MaxFrame)
                                
                                ElementX_j = this.A * XgYg(:,j,1);
                                ElementY_j = this.A * XgYg(:,j,2);
                                
                                Fx = DisplacementFields{i,1}{j};
                                Fy = DisplacementFields{i,2}{j};
                                
                                for e = 1:size(this.Elements, 1)
                                        Ux_i = ElementX_i(this.NumElementParameters*(e - 1) + 1:this.NumElementParameters*(e - 1) + this.NumElementParameters);
                                        Uy_i = ElementY_i(this.NumElementParameters*(e - 1) + 1:this.NumElementParameters*(e - 1) + this.NumElementParameters);
                                        Ux_j = ElementX_j(this.NumElementParameters*(e - 1) + 1:this.NumElementParameters*(e - 1) + this.NumElementParameters);
                                        Uy_j = ElementY_j(this.NumElementParameters*(e - 1) + 1:this.NumElementParameters*(e - 1) + this.NumElementParameters);
                                        
                                        Points_i = this.CalculatePoint(MaterialX(:), MaterialY(:), Ux_i, Uy_i);
                                        Points_j = this.CalculatePoint(MaterialX(:), MaterialY(:), Ux_j, Uy_j);

                                        Xj_ = Points_i(:,1) + feval(Fx, [Points_i(:,1), Points_i(:,2)]);
                                        Yj_ = Points_i(:,2) + feval(Fy, [Points_i(:,1), Points_i(:,2)]);
                                        
                                        weight = 1 / DisplacementFields{i,4}{j};

                                        metric = metric + (hypot(Points_j(:,1) - Xj_, Points_j(:,2) - Yj_).^2 .* weight);
                                end
                        end
                end
        end
end