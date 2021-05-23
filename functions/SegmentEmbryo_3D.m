function [toThreshold,RFP_FTL_RGB,Stats_table] = SegmentEmbryo_3D(toThreshold,LoGradius, Function, parameters,show,varargin)
    %RFP_FTL = zeros(size(toThreshold));
    RFP_FTL_RGB = zeros(size(toThreshold,1),size(toThreshold,2),3,size(toThreshold,3));
    %Stats = cell(size(toThreshold,3),1);
    Stats_table = cell(size(toThreshold,4),1);
    for f=1:size(toThreshold,4)
        disp(['segmenting f',num2str(f),'...'])
        toThres = FiltGlobalNorm(toThreshold(:,:,:,f),LoGradius);
        %[T_L T_L_RGB Stats] = FiltThresLab(img,LoGradius, level, areaopen)
        %[T_L T_L_RGB Stats] = ThresLabNBs(toThreshold,remove1, diskSize, WatershedParameter,remove2)
        [toThreshold(:,:,:,f),Stats_table{f}] = Function(toThres,parameters,varargin);
        %[toThreshold(:,:,:,f),Stats_table{f}] = ThresLabNBs_3D(toThres,num2cell([RemoveSize DiskSize WatershedParameter RemoveSize ThresLevel Thicken]),[XYRes,XYRes,ZRes]);
        %figure;subplot(131); imagesc(max(toThreshold,[],3))
                %subplot(132); imagesc(max(toThres,[],3))

    end
    RFP_FTL_max = MAX_proj_3D(toThreshold);
    for f=1:size(RFP_FTL_max,3)
        RFP_FTL_RGB(:,:,:,f)=label2rgb(RFP_FTL_max(:,:,f),'jet', 'k', 'shuffle');
    end
    
    if strcmp(show,'on')==1;
        D=zeros(size(RFP_FTL_max,1),size(RFP_FTL_max,2),1,size(RFP_FTL_max,3));
        D(:,:,1,:)=RFP_FTL_max(:,:,:);
        montage(D, [0 1]);
        D(:,:,1,:)=RFP_FTL_max(:,:,:)+1;
        mov = immovie(RFP_FTL_RGB);
        implay(mov)
    end
end