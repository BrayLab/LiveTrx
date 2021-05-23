function [F] = FiltGlobalNorm(toThreshold,LoGradius)
    F = zeros(size(toThreshold,1),size(toThreshold,2),size(toThreshold,3));
    h=-fspecial('log',round(10*1),LoGradius);
    for f=1:size(toThreshold,3)
        F(:,:,f) = fourierFilterWithSymmetricBoundaryConditions(toThreshold(:,:,f),h);
    end
    F(:,:,:)=F(:,:,:)./max(max(max(F(:,:,:))));
    disp('filtered ...')
    
end