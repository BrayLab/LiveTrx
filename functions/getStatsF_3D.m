function [Stats_GFP] = getStatsF_3D(varargin)
    RFP_FTL_tracked = varargin{1};
    F_max = varargin{2};
    
    Stats_GFP = cell(size(RFP_FTL_tracked,4),1);
    for f = 1:size(RFP_FTL_tracked,4)
        Stats_GFP{f,1} = regionprops('table',RFP_FTL_tracked(:,:,:,f),F_max(:,:,:,f),'PixelIdxList','PixelList','Centroid','MaxIntensity','MeanIntensity','PixelValues','Area');
    end
end