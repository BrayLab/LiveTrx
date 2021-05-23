function[TL,RFP_FTL_RGB,Stats] = ThresNucMemb_3D(toThreshold, parameters, refFilt)
%%
[GaussFilt, ThresLevel, Width, Height, Slices, XYRes, ZRes, Bits, WatershedParameter, Thicken, RemoveSize, RemoveSizeHigh] = parameters{:};
tic
% Thres = 0.05
% GaussFilt = 1;
%PropEmpty = 0.01
% WatershedParameter = 0.25;
% Thicken = 0;
% RemoveSize = 10; 
% RemoveSizeHigh = 200; %in um^3

%ZRes = 1
    %toThreshold = ReadSingleStack(reader,Channels,Slices,Frames0,Flip,From,ChannelToTrack,f);
%toThreshold = imresize3(toThreshold,[Width,Height,round(Slices*ZRes/XYRes)]);
%ZResOld = ZRes; ZRes = XYRes;
 
%Filt = Filter_3D(toThreshold, MedFilt, 'off');
Filt = imgaussfilt3(toThreshold,GaussFilt)./(2^Bits-1);
% adjust contrast based on mean and SD (no extra) to make initial
% segmentation. then active contours seeded by this wo applying increase in
% contrast to input of AC
FiltC = imhistmatchn(Filt,refFilt);
% Mean = mean2(max(Filt,[],3));
% SD = std2(max(Filt,[],3));
% %Mean = mean(toThreshold(:))
% %SD = std(toThreshold(:))
% Low = Mean-(SD*3); High = Mean+(SD*2);
% if Low < 0; Low = 0; end
% if High > 1; High = 1; end
% FiltC = imadjustn(Filt,[Low, High],[0; 1]);
%Thres = graythresh(FiltC)
BW = FiltC > ThresLevel; %0.05
%BW = Filt > Thres; %0.05


% BW = zeros(size(toThreshold));
% for z = 1:size(toThreshold, 3)
%     BW(:,:,z) = toThreshold(:,:,z) > 0.1*255;
% end

AC = activecontour(Filt, BW, 50);
ACi = ~AC;
%GDW = graydiffweight(toThreshold,max(toThreshold(:))-mean(toThreshold(:)));

%if 1==0
%slices with mostly 1 pixels (above and below nuc because inverted)
        %will be changed to 0. slices from vit memb to firstnuc also
        %changed to 0
        P = permute(ACi,[3,2,1]);
        P = sum(P,3); P = 1-sum(P,2)./max(sum(P,2));
        % if at least one of the first 3 slices is vit membrane from first
        % to first nuc replaced by ones (0 in matrix)
        VitMemb = P >= 0.8;
        FirstNuc = find(islocalmin(P) & (~VitMemb));
        if sum(VitMemb(1:10)) >= 1 & ~isempty(FirstNuc)% if at least one of first 10 slices is vitmemb
            P(1:FirstNuc(1)) = 1;
        end
        plot(P); hold on
        plot(islocalmin(P) & (~VitMemb))
        drawnow
 %replace above and below nuc, and vit memb if there with 1s before wat       
 ACi(:,:,P==1) = 1;
%end        

        %AC = flip(AC,3);


          if isnan(WatershedParameter)~=1
            try
                %Voxel = varargin{1};
                Voxel = [XYRes,XYRes,ZRes];
                D = bwdistsc(~ACi,Voxel);
            catch
                D = bwdist(~ACi);
            end
            D = -D;
            %D(~L) = -Inf;
            D2 = imhmin(D,WatershedParameter);
            %D2=D;
            D2(~ACi) = -Inf;
            L2 = watershed(D2, 18);
            L3 = bwareaopen(L2-1,round(RemoveSize*(XYRes)^3),6);
        else
            L3 = ACi;
        end
        
            TL = bwlabeln(L3,6);
  
        
            %TL = flip(TL,3);

        
        cmap = jet(max(TL(:))+1);
        cmap(1,:) = [0,0,0];
        %imagesc(TL(:,:,60)); colormap(cmap)
        

        Stats_table = regionprops('table',TL,'Area','Centroid',...
            'BoundingBox');
        
        %Stats_table = regionprops3('table',TL,'Volume','Centroid',...
        %    'Orientation','SubarrayIdx','PrincipalAxisLength','SurfaceArea');
        
        %plot(ones(size(Stats_table.Area))+rand(size(Stats_table.Area)),Stats_table.Area,'.')

        %histogram(Stats_table.Area.*(XYRes^3), 'BinWidth',10)
                %histogram(Stats_table.BoundingBox(1,:), 'BinWidth',10)

        Stats_table.Label = [1:length(Stats_table.Area)]';
        ToRemove = find((Stats_table.Area.*(XYRes^3) < RemoveSize) | (Stats_table.Area.*(XYRes^3) > RemoveSizeHigh) | (Stats_table.BoundingBox(:,1) < 1));
        Stats_table(ToRemove,:) = [];
        PixelsToRemove = ismember(TL,ToRemove);
        %REmoved = zeros(size(TL)); REmoved(PixelsToRemove) = 1;
        TL(PixelsToRemove) = 0;
        close all
        %SE = strel('sphere',3);
        %TL = imdilate(TL,SE);
        if Thicken~=0
            for z=1:size(TL,3)
                TL(:,:,z) = bwmorph(TL(:,:,z),'thicken',Thicken);
            end
            TL = bwlabeln(TL,6);
        end



        Stats_table = regionprops('table',TL,'Area','Centroid',...
            'BoundingBox');
        
        % second round of filtering by size because dilation might have
        % created new objects (FIND A BETTER WAY!!!)
     if ~isempty(Stats_table)
        Stats_table.Label = [1:length(Stats_table.Area)]';
        ToRemove = find((Stats_table.Area.*(XYRes^3) < RemoveSize) | (Stats_table.Area.*(XYRes^3) > RemoveSizeHigh) | (Stats_table.BoundingBox(:,1) < 1));
        Stats_table(ToRemove,:) = [];
        PixelsToRemove = ismember(TL,ToRemove);
        %REmoved = zeros(size(TL)); REmoved(PixelsToRemove) = 1;
        TL(PixelsToRemove) = 0;
        close all
     end
        TL = bwlabeln(TL,6);
        
        %Overlay = Filt + double(TL~=0);

        imagesc(max(TL,[],3)); colormap(cmap);
        drawnow

        Stats_table = regionprops('table',TL,'Area','Centroid',...
            'BoundingBox','PixelIdxList');
        
        Stats{1,1} = Stats_table;
        %Stats_table = regionprops3(TL,'Volume','Centroid',...
        %    'Orientation','SubarrayIdx','PrincipalAxisLength','SurfaceArea');
        
        RFP_FTL_RGB = label2rgb(max(TL,[],3),'jet', 'k', 'shuffle');
        
        %histogram(Stats_table.PrincipalAxisLength(:,1)); hold on
        %histogram(Stats_table.PrincipalAxisLength(:,2))
        %histogram(Stats_table.PrincipalAxisLength(:,3))


toc        
end
        %%
        
        
        
        