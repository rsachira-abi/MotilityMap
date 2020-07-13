
function ssd = NonlinearMinimization (this, PxPy, GivenPx, GivenPy, GivenRange, CurrentRange, ConvertionMatrix)
    PxPy = reshape(PxPy, [], length(CurrentRange), 2);
    
    Px = cat(2, GivenPx, PxPy(:,:,1));
    Py = cat(2, GivenPy, PxPy(:,:,2));
    
    TotalRange = [GivenRange, CurrentRange];
    
    ssd = 0;
    
    EulerPointsX = ConvertionMatrix * Px;
    EulerPointsY = ConvertionMatrix * Py;
    
    for i = 1:size(EulerPointsX, 2)
        
        I = TotalRange(i);
        if (I > size(this.DisplacementFields, 1))
            break;
        end
        
        [Lia, Locb] = ismember(CurrentRange, this.DisplacementFields{I,3});
        
        for k = find(Lia)
            j = Locb(k);
            J = k + length(GivenRange);
            
            Fx = this.DisplacementFields{I,1}{j};
            Fy = this.DisplacementFields{I,2}{j};
            
            Gx = Fx(EulerPointsX(:,i), EulerPointsY(:,i));
            Gy = Fy(EulerPointsX(:,i), EulerPointsY(:,i));
            
            weight = this.DisplacementFields{I,4}(j) / this.MaxWeight;
            %weight = 1 / j;
            
            error = (EulerPointsX(:,i) + Gx - EulerPointsX(:,J)).^2 + (EulerPointsY(:,i) + Gy - EulerPointsY(:,J)).^2;
            
            ssd = ssd + (sum(error) .* weight);
        end
    end
end



