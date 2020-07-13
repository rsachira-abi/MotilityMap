% Helper functions for unit tests
classdef Utils
    
    methods(Static, Access = public)
        
        function [template, target] = generateShiftedPair (shiftX, shiftY)
            I = imread('noise.png');
            I = im2double(rgb2gray(I));

            [baseMeshX, baseMeshY] = meshgrid(100:400, 100:400);
            
            shiftedMeshX = baseMeshX - shiftX;
            shiftedMeshY = baseMeshY - shiftY;
            
            template = interp2(I, baseMeshX, baseMeshY, 'spline');
            target = interp2(I, shiftedMeshX, shiftedMeshY, 'spline');
        end
        
    end
end