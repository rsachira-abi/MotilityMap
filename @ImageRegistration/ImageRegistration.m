%
% c - indicate variables in residing in the gpu (cuda)

classdef ImageRegistration < handle
    properties (Constant)
        PeakThreshold = 0.75;
        GradKernel = [22/252,-67/252,-58/252,0,58/252,67/252,-22/252]; % SGD
    end
    
    properties (Access = private)
        %% Configuration
        HalfCCSize = int32(64); % int32(254); % int32(128); % int32(64);
        CCSize = 129; %509; % 257; % 129;
        SquareSize = int32(12);
        
        %% Other variables
        cHamming2D = [];
        cHann2D = [];
        
        Points = [];
        cPaddedPoints = [];
        cPaddedTemplateGrad = [];
        cNextPaddedTemplate = [];
        
        A = [];
        cA = [];
    end
    
    %% Public methods
    methods (Access = public)
        
        function this = ImageRegistration (template, points, HalfCCSize)
			if (nargin == 3)
				this.HalfCCSize = int32(HalfCCSize);
				this.CCSize = double(this.HalfCCSize * 2 + 1);
				this.SquareSize = int32(sqrt(double(this.CCSize / 4)));
			end
		
            hamming2D = hamming(this.CCSize) * hamming(this.CCSize)';
            this.cHamming2D = hamming2D; %gpuArray(hamming2D);
            
            hann2D = hann(this.CCSize) * hann(this.CCSize)';
            this.cHann2D = hann2D; %gpuArray(hann2D);
            
            % Matrix for plane fitting.
            %len = 61;
            len = this.CCSize - 6 + 1 - (2 * round(this.CCSize * 0.25)) + 1;
            [mX, mY] = meshgrid(1:len, 1:len);
            this.A = cat(2, mX(:), mY(:), ones(size(mX(:))));
            this.cA = this.A; %gpuArray(this.A);
            
            if (nargin > 0)
                templateGrad = this.grad(template);
                [PaddedTemplateGrad, padded_points] = this.padArray(templateGrad, points);
                
                this.Points = points;
                this.cPaddedPoints = padded_points; %gpuArray(padded_points);
                this.cPaddedTemplateGrad = PaddedTemplateGrad; %gpuArray(PaddedTemplateGrad);
            end
        end
        
        [DispX, DispY] = registerImage (TemplateGrad, TargetGrad, Points);
        
        grad_img = grad (this, image);
        
        function [pad_image, pad_points] = padArray(this, image, points)
            pad_image = padarray(image, [this.CCSize, this.CCSize]);
            pad_points = points + this.CCSize;
        end
        
        function sz = getTemplateSize (this)
            sz = size(this.cPaddedTemplateGrad);
        end
    end
    
    %% Private methods
    methods (Access = public)
        
        [fieldsX, fieldsY, error] = calcDisplacementFields (this, target_grads)
        
        [intX, intY, intError] = calcIntegerShift (this, cPaddedTemplateGrad, cPaddedTargetGrad, cPoints, cDispX, cDispY)
        
        [subX, subY] = calcSubpixelShift (this, cPaddedTemplateGrad, cPaddedTargetGrad, cPoints, cDispX, cDispY)
        
    end
    
    %% CUDA functions
    methods (Static = true, Access = public)
        
        I = ImportRaw(filename);
    end
end