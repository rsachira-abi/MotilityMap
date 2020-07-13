
function [E, F] = GenerateConversionMatrices (this, MaterialPoints)
    
    for x = 0:this.n
        for y = 0:this.n
            r = x * (this.n + 1) + y + 1;

            for i = 0:this.n
                for j = 0:this.n
                    c = i * (this.n + 1) + j + 1;
                    
                    E(r,c) = sum(this.N_(i, this.k, MaterialPoints(:,1)) .* this.N_(j, this.k, MaterialPoints(:,2))...
                        .* this.N_(x, this.k, MaterialPoints(:,1)) .* this.N_(y, this.k, MaterialPoints(:,2)));
                end
            end
            
            F(r,:) = (this.N_(x, this.k, MaterialPoints(:,1)) .* this.N_(y, this.k, MaterialPoints(:,2)))';
        end
    end
end