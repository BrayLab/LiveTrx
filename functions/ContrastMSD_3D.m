function [RFPtoThreshold] = ContrastMSD_3D(RFPtoThreshold, InputLow, InputHigh,Bits,show)
    for f=1:size(RFPtoThreshold,4)
        disp(['contrasting f',num2str(f),'...'])
        toThreshold = RFPtoThreshold(:,:,:,f)./(2^Bits-1);
        Mean = mean2(MAX_proj_3D(toThreshold));
        SD = std2(MAX_proj_3D(toThreshold));
        %Mean = mean(toThreshold(:))
        %SD = std(toThreshold(:))
        try
            Low = Mean-SD; High = Mean+SD;
            if Low < 0; Low = 0; end
            if High > 1; High = 1; end
            %for z=1:size(toThreshold,3)
                % first adjust for Mean and SD and then for extra contrast
                % if specified as argument (most times InputLow = 0 and
                % InputHigh = 1
             %   toThreshold(:,:,z) = imadjust(toThreshold(:,:,z),[Low, High],[0; 1]);
              %  toThreshold(:,:,z) = imadjust(toThreshold(:,:,z),[InputLow,InputHigh],[0; 1]);
            %end
            toThreshold = imadjustn(toThreshold,[Low, High],[0; 1]);
            toThreshold = imadjustn(toThreshold,[InputLow,InputHigh],[0; 1]);
         
        end
        RFPtoThreshold(:,:,:,f) = toThreshold;
    end

    if strcmp(show,'on')==1;
        D=zeros(size(RFPtoThreshold,1),size(RFPtoThreshold,2),1,size(RFPtoThreshold,4));
        %D(:,:,1,:)=RFPtoThreshold(:,:,1,:);
        %montage(D, [0 1]);
        D(:,:,1,:)=RFPtoThreshold(:,:,floor(size(RFPtoThreshold,3)/2),:);
        montage(D, [0 1]);
        %D(:,:,1,:)=RFPtoThreshold(:,:,size(RFPtoThreshold,3),:);
        %montage(D, [0 1]);
    end
end