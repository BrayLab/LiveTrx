function[] = mainSegTrackMS2_3D_fromBF_combined(Parameters,varargin)
%
Table2Vars(Parameters);
Flip = str2double(strsplit(Flip,','));

mkdir([Path,File,Name]) %uncomment
PathToSave = [Path,File,Name,File]; %uncomment
show = 'off';
reader = bfGetReader([Path,File]);
[Bits,Width,Height, Channels, Slices, Frames, Frames0, XRes, YRes, ZRes,Zoom,TimeRes,Settings] = readMetadataBFOME(reader);
try
    XYRes=mean(round(mean([XRes,YRes]),2));
    ZRes = mean(round(ZRes,2));
end
if strcmp(To,'NA')==1 | To == inf;  To=Frames-1;
elseif isnumeric(To) == 0; To = str2num(To); end;
Frames = To-From+1;

if strcmp(Channel1,'nucmemb') | strcmp(Channel2,'nucmemb')
    NucMemb = 1;
else
    NucMemb = 0
end
if strcmp(Channel1,'memb') | strcmp(Channel2,'memb')
    Membranes = 1;
else
    Membranes = 0
end
%% get settings from arguments or use default
InputLow = 0; InputHigh = 1;
try 
    Settings = varargin{1};
    cellfun(@(x) evalin('caller',x), Settings)
    DiskSize = round(1.5./XYRes); %in px, not used? but error if it doesnt exist
    fileID = fopen([PathToSave,'_Settings_SegTrack.txt'],'a')
    fprintf(fileID, '\nupdated on \t%s . used settings:\n',date)
    fprintf(fileID,'%s \n',Settings{:})
    fclose(fileID);
    SplitShortLong = 60; % ignore
    disp('loaded settings')
catch
    if NucMemb || Membranes
        Divisions = 1;
        Spots = 0;
        ThresLevel = 0.05 %FIXED THRESHOLD LEVEL FOR INITIAL SEGMENTATION
        GaussFilt = 1; %GAUSSIAN FILTERING RADIUS
        WatershedParameter = 0.25; % AMOUNT OF WATERSHEDING (when lower more watershed, when higher less)
        Thicken = 0; % leave at 0
        RemoveSize = 10; % objects smaller than this (in um^3) will be deleted
        RemoveSizeHigh = 150; % objects bigger than this (in um^3) will be deleted
        %Divisions = 1; %leave at 1. how to fix conflicts when two nuclei match to 1. it will add new labels for both of them but most likely match can be fixed afterwards
        Distance = 4; % maximum distance (in um) for finding a closest neighbour
        MaxN = 5; % maximum # of frames back that tries to find a closest neighbour, if no one found, a new label is created
        PrintInfo = 0;
    else
        Divisions = 1;
        Spots = 0;
        MedFilt = 3; %%3 in px
        GaussFilt = 0; %GAUSSIAN FILTERING RADIUS
        InputLow = 0; %%0.1
        InputHigh = 1; % was 1
        LogRadius = 4./XYRes*0.4; %0.4 from HG's scripts
        ThresLevel = 0.05;
        DiskSize = round(1.5./XYRes); %in px, not used?
        WatershedParameter = 1; % was 1
        Thicken = round(1./XYRes); %in px
        RemoveSize = round(10./(XYRes^2.*ZRes)); %10 um^3 in px volume. nc14 nuc is aprox 33 um^3
        Distance = 6; %nuclei allowed to move 6um
        MaxN = 2;
        PrintInfo = 0;
    end
end

%% INITIALIZE

Segmented_RGB = zeros(Height, Width, 3,Frames);
Stats_tracked = cell(Frames,1);
MaxL = 1; 



ChannelToTrack = 0; %TO DO, when two transcription channels
if strcmp(Channel1,'his') | strcmp(Channel1,'nucmemb') | strcmp(Channel1,'memb')
    ChannelToTrack = 1;
end
if strcmp(Channel2,'his') | strcmp(Channel2,'nucmemb') | strcmp(Channel2,'memb')
    ChannelToTrack = 2;
end
if strcmp(Channel3,'his') | strcmp(Channel3,'nucmemb') | strcmp(Channel3,'memb')
    ChannelToTrack = 3;
end
%ChannelToTrack = 1; % ###############################################

ChannelTrx = [0,0,0]; %TO DO, when two transcription channels. DONE 28/02/2021
if strcmp(Channel1,'MCP') || strcmp(Channel1,'PCP')
    ChannelTrx(1) = 1;
end
if strcmp(Channel2,'MCP') || strcmp(Channel2,'PCP')
    ChannelTrx(2) = 1;
end
if strcmp(Channel3,'MCP') || strcmp(Channel3,'PCP')
    ChannelTrx(3) = 1;
end

ChannelTrx = find(ChannelTrx);

ChannelMean = [0,0,0];
if strcmp(Channel1,'mean')
    ChannelMean(1) = 1;
    if ChannelToTrack == 0; ChannelToTrack = 1; end
end
if strcmp(Channel2,'mean')
    ChannelMean(2) = 1;
end
if strcmp(Channel3,'mean')
    ChannelMean(3) = 1;
end
ChannelMean = find(ChannelMean);

ChannelsDescription = {Channel1,Channel2,Channel3};

if NucMemb || Membranes
    toTrack = zeros(Height, Width, round(Slices*ZRes/XYRes),Frames);
    disp(['filtering ref frame, f',num2str(nc14),'...']);
    toThreshold = ReadSingleStack(reader,Channels,Slices,Frames0,Flip,From,ChannelToTrack,nc14);
    toThreshold = imresize3(toThreshold,[Width,Height,round(Slices*ZRes/XYRes)]);
    refFilt = imgaussfilt3(toThreshold,GaussFilt)./(2^Bits-1);
else
    toTrack = zeros(Height, Width, Slices,Frames);
end

%change for when both channels have transcription

%% SEGMENTATION

for f = 1:Frames
    disp(['reading and segmenting f',num2str(f),'...']);
    toThreshold = ReadSingleStack(reader,Channels,Slices,Frames0,Flip,From,ChannelToTrack,f);
    if Membranes
        toThreshold = imresize3(toThreshold,[Width,Height,round(Slices*ZRes/XYRes)]);
    end
    if NucMemb || Membranes
        % resize to make isotropic voxels
        toThreshold = imresize3(toThreshold,[Width,Height,round(Slices*ZRes/XYRes)]);
        parameters = num2cell([GaussFilt, ThresLevel, Width, Height, Slices, XYRes, XYRes, Bits, WatershedParameter, Thicken, RemoveSize, RemoveSizeHigh]);
        [Dummy1,Dummy2,Dummy3] = ThresNucMemb_3D(toThreshold, parameters, refFilt);
    else
        %added filter to remove vitelline membrane when signal is very low
        if strcmp(Channel1,'mean') == 1 && ChannelToTrack ==1
            toThreshold = imgaussfilt3(toThreshold,1);
            InputHigh = 0.5;
%             MaskVitMemb = imbinarize(toThreshold./max(toThreshold(:)));
%             SE =strel('sphere',2);
%             MaskVitMemb = imdilate(MaskVitMemb,SE);
%             toThreshold(MaskVitMemb) = 0;
        end
        [toThreshold] = Filter_3D(toThreshold, MedFilt, 'off');
        if GaussFilt ~= 0
        	toThreshold = imgaussfilt3(toThreshold,GaussFilt);
        end
        [toThreshold] = ContrastMSD_3D(toThreshold,InputLow, InputHigh,Bits,'off');
        [Dummy1 ,Dummy2 , Dummy3] = SegmentEmbryo_3D(toThreshold,LogRadius, @ThresLabNBs_3D, num2cell([RemoveSize DiskSize WatershedParameter RemoveSize ThresLevel Thicken]),'off',[XYRes,XYRes,ZRes]);
    end
    toTrack(:,:,:,f) = Dummy1;
    Segmented_RGB(:,:,:,f) = Dummy2(:,:,:,1);
    Stats_tracked{f,1} = Dummy3{1,1}; 
    MaxL = max(MaxL, max(max(max(Dummy1))));
end

% new ZRes after isotropic voxels, old one stored
if NucMemb || Membranes
    ZResOld = ZRes; ZRes = XYRes;   
end


%% TRACKING
% cant put into function because it would duplicate toTrack when modifying
% it and it crashes
    cmap = jet(1000000);
    cmap_shuffled = cmap(randperm(size(cmap,1)),:);
    newLabel=MaxL+1;
    Stats_tracked{1,:}.Label = (1:size(Stats_tracked{1,:},1))';
    if Divisions
            Stats_tracked{1,:}.Parent = zeros(size(Stats_tracked{1,:},1),1);
            Stats_tracked{1,:}.Daughters = zeros(size(Stats_tracked{1,:},1),2);
    end
    toReplace = Stats_tracked{1,:}(1,:); 
    for r = 1:size(toReplace,2); try; toReplace(1,r) = table(NaN);end; end
     
    TrackedBoundL = zeros(Height, Width,Frames);
    TrackedBoundL(:,:,1) = MAX_proj_3D(toTrack(:,:,:,1));
    for f=2:size(toTrack,4)
        disp(['frame ',num2str(f-1),' to ',num2str(f)])
        F1_0 = toTrack(:,:,:,f);
        F1_t = zeros(size(toTrack(:,:,:,f)));
        FromT = f-MaxN; if FromT < 1; FromT = 1; end
        SubStats = Stats_tracked(FromT:f,1);
        [F1_t, ~, SubStats, newLabel] = TrackWithDivisions_3D(SubStats,F1_0,F1_t, cmap_shuffled,f,Distance,newLabel,MaxN,Divisions,toReplace,XYRes, ZRes);  %for others   
        toTrack(:,:,:,f) = F1_t(:,:,:);
        Stats_tracked(FromT:f,1) = SubStats;
        TrackedBoundL(:,:,f) = MAX_proj_3D(F1_t);
    end
    %[toTrack , ~, Stats_tracked] = Tracking_3D(toTrack, Stats_table, cmap_shuffled, Distance,MaxN,Divisions,'off','off');
    %clear Stats_table
    [BoundariesBW, Boundaries_RGB, ~] = BoundariesTracked_3D(toTrack,cmap_shuffled,show);
    %[Boundaries_RGB] = double(TimeStamp(uint8(Boundaries_RGB),TimeRes,nc14,Delay));
    WriteRGB(Boundaries_RGB, PathToSave, '_segmented_tracked_boundaries_RGB.tiff','none')
    %delete empty lines (needed for tracking) before saving
    for f= 1:size(Stats_tracked,1);
        if ~isempty(Stats_tracked{f,1})
            Stats_tracked{f,1}(isnan(Stats_tracked{f,1}.Label)|Stats_tracked{f,1}.Area==0,:) = [];
        end
    end
    save([Path,File,Name,File,'_Stats_tracked.mat'],'Stats_tracked','-v7.3');

%% clean objects tracked for less than minimum after tracking and deleting 
% them from toTrack (CHANGE SO THAT IT DOES IT AFTER FIX TRACKING)

if NucMemb || Membranes
    minNumb = 5
    [AllF] = MergeAll(Stats_tracked, TimeRes);
    AllF = splitvars(AllF);
        %
    Labels = unique(AllF.Label); LabelsOld = Labels;
    [Area] = Reshape(AllF,Frames,Labels,'Area','Label');
        
    noNaNnum = sum(~isnan(Area),1);
    ToDelete = (noNaNnum <= minNumb);
    % delete tracks less than minNumb
    LabelsToDelete = Labels(ToDelete);
    Labels(ToDelete) = [];
    Area(:,ToDelete) = [];

    for f = 1:size(toTrack,4)
        disp(num2str(f))
        TL = toTrack(:,:,:,f);
        PixelsToRemove = ismember(TL,LabelsToDelete);
        TL(PixelsToRemove) = 0;
        toTrack(:,:,:,f) = TL;
    end
    [~, Boundaries_RGB_clean, ~] = BoundariesTracked_3D(toTrack,cmap_shuffled,show);
    %[Boundaries_RGB] = double(TimeStamp(uint8(Boundaries_RGB),TimeRes,nc14,Delay));
    WriteRGB(Boundaries_RGB_clean, PathToSave, '_segmented_tracked_boundaries_RGB_clean.tiff','none')
    %printLabels_new(Boundaries_RGB,Stats_table,2,'off', PathToSave, '_segmented_tracked_info.tiff','none')
end
 %%   
if Divisions
    [~,FTL_tracked_divisions_RGB,Phylogeny, Families] = PhylogenyDivisions(Stats_tracked,TrackedBoundL,1, Frames,cmap_shuffled,'on');
    writetable(cell2table(Phylogeny),[Path,File,Name,File,'_phylogeny.txt']);
    writetable(Families,[Path,File,Name,File,'_families.txt']);
    clear Phylogeny; clear Families;
    %[FTL_tracked_divisions_RGB] = double(TimeStamp(uint8(FTL_tracked_divisions_RGB),TimeRes,nc14,Delay));
    %WriteRGB(FTL_tracked_divisions_RGB, PathToSave, '_segmented_tracked_divisions_RGB.tiff','none')
    clear FTL_tracked_divisions_RGB;
end
boundariesL = BoundariesBW .* TrackedBoundL;
%Write16b(boundariesL, PathToSave, '_boundariesL.tiff', 'none');
%clear BoundariesBW; clear TrackedBoundL;
    %
   
%% rescue toTrack from Stats_GFP
%         load([PathToSave, '_Stats_tracked.mat']);
% toTrack = zeros(Height,Width, Slices, Frames);
%     for f = 1:length(Stats_tracked)
%         disp(num2str(f))
%          SingleFrame = Stats_tracked{f};
%          Tracked3D = zeros(Height,Width,  Slices);
%          for n = 1:size(SingleFrame,1)
%              Tracked3D(SingleFrame.PixelIdxList{n}) =  SingleFrame.Label(n);
%          end
%          toTrack(:,:,:,f) = Tracked3D;
%     end
   %% MEASURE FLUORESCENCE IN MCP/PCP/mean channels
ChannelsMeasure = [ChannelTrx,ChannelMean];
if ~isempty(ChannelsMeasure)
for CTrx = 1:length(ChannelsMeasure)
    MAX_G = zeros(Height, Width,Frames);
    MAX_replacedF = zeros(Height, Width,Frames);
    Stats_Trx = cell(Frames,1);
    StatsSpots = cell(Frames,1);
    if Spots ~= 0
        MAX_GSeg = zeros(Height, Width,Frames);
        MAX_GGaussF = zeros(Height, Width,Frames);
        MAX_meanFGauss = zeros(Height, Width,Frames);
    end
    maxWorkers = 6;
    if Spots ~= 0
        try
            parpool(maxWorkers);  % 6 is the number of cores the Garcia lab server can reasonably handle per user at present.
        catch
            try
                parpool; %in case there aren't enough cores on the computer 
            catch
            end
            %parpool throws an error if there's a pool already running. 
        end
    end
    for f = 1:Frames
        disp(['reading f',num2str(f),'...']);
        G = ReadSingleStack(reader,Channels,Slices,Frames0,Flip,From,ChannelsMeasure(CTrx),f);
        [Stats] = getStatsF_3D(toTrack(:,:,:,f), G); 
        [Stats] = PrintF_3D(Stats,Path,File,Name,0);
        MAX_G(:,:,f) = MAX_proj_3D(G);
        % only run spot tracking when specified and when tracking
        % transcription
        if Spots ~= 0 & ~strcmp(ChannelsDescription{ChannelsMeasure(CTrx)},'mean')
            GSeg = zeros(Height, Width,Slices);
            GGaussF = zeros(Height, Width,Slices);
            Threshold = 150;
            Shadows = 0;
            [SSpots] = segmentSpots_JF(G, Threshold,Shadows, XYRes);
            %SpotsAll{f} = SSpots;
            if numel(fieldnames(SSpots)) ~= 0
                Positions = [[SSpots.Fits.CentX]', [SSpots.Fits.CentY]', [SSpots.Fits.brightestZ]'];
                Intensities = [SSpots.Fits.MaxInt]';
                Indices = [sub2ind(size(G), Positions(:,2),Positions(:,1),Positions(:,3))];
                GSeg(Indices) = 1;
                GGaussF(Indices) = Intensities;
            end
            SpotsLabelled = bwlabeln(GSeg,6);
            [StatsS] = getStatsF_3D(SpotsLabelled, GGaussF); 
            StatsSpots{f} = StatsS{1,1};

            GSeg = max(GSeg,[],3);
            GSeg = bwmorph(GSeg, 'thicken',3);
            GSeg = GSeg*4095;
            MAX_GSeg(:,:,f) = GSeg;
            MAX_GGaussF(:,:,f) = max(GGaussF,[],3);

            [StatsGauss] = getStatsF_3D(toTrack(:,:,:,f), GGaussF); 
            [StatsGauss] = PrintF_3D(StatsGauss,Path,File,Name,Spots);
            try
                Stats{1,1}.MaxGauss = StatsGauss{1,1}.Max;
                Stats{1,1}.PositionGaussX = StatsGauss{1,1}.SpotPositionX;
                Stats{1,1}.PositionGaussY = StatsGauss{1,1}.SpotPositionY;
                Stats{1,1}.PositionGaussZ = StatsGauss{1,1}.SpotPositionZ;
                [ReplacedSpots,Stats] = ReplaceLabelsbyF_3D(toTrack(:,:,:,f), Stats, 0,2^(Bits+1)-1,'MaxGauss');
                MAX_meanFGauss(:,:,f) = max(ReplacedSpots,[],3);
                Stats{1,1}.MaxGauss2 = StatsGauss{1,1}.Max2;
                Stats{1,1}.Position2GaussX = StatsGauss{1,1}.Spot2PositionX;
                Stats{1,1}.Position2GaussY = StatsGauss{1,1}.Spot2PositionY;
                Stats{1,1}.Position2GaussZ = StatsGauss{1,1}.Spot2PositionZ;
            end
        end 
         % ############################################### replace by max
         % or mean intensity depending on what's being tracked
        if ~strcmp(ChannelsDescription{ChannelsMeasure(CTrx)},'mean')
            [toTrackFrameReplaced,Stats] = ReplaceLabelsbyF_3D(toTrack(:,:,:,f), Stats, 0,2^Bits-1,'Max');
        else
            [toTrackFrameReplaced,Stats] = ReplaceLabelsbyF_3D(toTrack(:,:,:,f), Stats, 0,2^Bits-1,'MeanIntensity');
        end
        Stats_Trx{f,1} = Stats{1,1};
        MAX_replacedF(:,:,f) = max(toTrackFrameReplaced,[],3);
    end
    clear G
%
%
    %clear toTrack
    Merged_replacedF_maxGFP = (MAX_G./(2.^(Bits-8))+ MAX_replacedF);
    Merged_replacedF_maxGFP = uint8(permute(Merged_replacedF_maxGFP,[1,2,4,3]));
    Merged_meanF_maxGFP_RGB = cat(3,Merged_replacedF_maxGFP,Merged_replacedF_maxGFP,Merged_replacedF_maxGFP);
    [Merged_meanF_maxGFP_RGB] = TimeStamp(Merged_meanF_maxGFP_RGB,TimeRes,nc14,Delay);
    Merged_replacedF_maxGFP(:,:,:) = Merged_meanF_maxGFP_RGB(:,:,1,:);
    Merged_replacedF_maxGFP = double(permute(Merged_replacedF_maxGFP,[1,2,4,3]));
    % file called maxF but will be mean when MeanIntensity is replaced
    Write8b(Merged_replacedF_maxGFP, PathToSave, ['_replacedF_maxGFP_',ChannelsDescription{ChannelsMeasure(CTrx)},'.tiff'])
    
        MAX_replacedF = uint8(permute(MAX_replacedF,[1,2,4,3]));
        Merged_meanF_maxGFP_RGB = cat(3,MAX_replacedF,MAX_replacedF,MAX_replacedF);
        [Merged_meanF_maxGFP_RGB] = TimeStamp(Merged_meanF_maxGFP_RGB,TimeRes,nc14,Delay);
        MAX_replacedF(:,:,:) = Merged_meanF_maxGFP_RGB(:,:,1,:);
        MAX_replacedF = double(permute(MAX_replacedF,[1,2,4,3]));
        Write8b(MAX_replacedF, PathToSave, ['_replacedF_',ChannelsDescription{ChannelsMeasure(CTrx)},'.tiff'])
    
    if Spots ~= 0
        Merged_replacedF_maxGFPGauss = (MAX_G./(2.^(Bits-8))+MAX_meanFGauss);
        Merged_replacedF_maxGFPGauss = uint8(permute(Merged_replacedF_maxGFPGauss,[1,2,4,3]));
        Merged_meanF_maxGFP_RGB = cat(3,Merged_replacedF_maxGFPGauss,Merged_replacedF_maxGFPGauss,Merged_replacedF_maxGFPGauss);
        [Merged_meanF_maxGFP_RGB] = TimeStamp(Merged_meanF_maxGFP_RGB,TimeRes,nc14,Delay);
        Merged_replacedF_maxGFPGauss(:,:,:) = Merged_meanF_maxGFP_RGB(:,:,1,:);
        Merged_replacedF_maxGFPGauss = double(permute(Merged_replacedF_maxGFPGauss,[1,2,4,3]));
        Write8b(Merged_replacedF_maxGFPGauss, PathToSave, ['_maxF_maxGFP_Gauss_',ChannelsDescription{ChannelsMeasure(CTrx)},'.tiff'])

        Merged_replacedF_maxGFPGauss = (MAX_G./(2.^(Bits-8))+MAX_GSeg./(2.^(Bits-8+1)));
        Merged_replacedF_maxGFPGauss = uint8(permute(Merged_replacedF_maxGFPGauss,[1,2,4,3]));
        Merged_meanF_maxGFP_RGB = cat(3,Merged_replacedF_maxGFPGauss,Merged_replacedF_maxGFPGauss,Merged_replacedF_maxGFPGauss);
        [Merged_meanF_maxGFP_RGB] = TimeStamp(Merged_meanF_maxGFP_RGB,TimeRes,nc14,Delay);
        Merged_replacedF_maxGFPGauss(:,:,:) = Merged_meanF_maxGFP_RGB(:,:,1,:);
        Merged_replacedF_maxGFPGauss = double(permute(Merged_replacedF_maxGFPGauss,[1,2,4,3]));
        Write8b(Merged_replacedF_maxGFPGauss, PathToSave, ['_maxF_maxGFP_Spots_',ChannelsDescription{ChannelsMeasure(CTrx)},'.tiff'])
    end

    if Spots == 2
        [Stats_Trx] = Retrack2Spots(Stats_Trx);
    end

% save used parameters
if NucMemb || Membranes
    parameters = table(Bits,Channels,Slices, sum(Frames0),Frames,Width,Height,From, To, GaussFilt,ThresLevel,WatershedParameter,Thicken,RemoveSize,RemoveSizeHigh,Distance,MaxN,TimeRes,XYRes,ZRes,Flip,Zoom);
    writetable(parameters,[Path,File,Name,File,'_parameters.txt']);
    Stats_Trx = Stats_tracked;
else
    parameters = table(Bits,Channels,Slices, sum(Frames0),Frames,Width,Height,From, To, MedFilt,InputLow,InputHigh,LogRadius,ThresLevel,WatershedParameter,DiskSize,Thicken,RemoveSize,Distance,MaxN,TimeRes,XYRes,ZRes,Flip,Zoom);
    writetable(parameters,[Path,File,Name,File,'_parameters.txt']);

disp('done')
save([Path,File,Name,File,'_Stats',ChannelsDescription{ChannelsMeasure(CTrx)},'.mat'],'Stats_Trx','-v7.3');

end
end
end
%clear toTrack
%% print "info" movie with label numbers in all cells
%% to print all boundaries
% 
if PrintInfo
    FTL_tracked_meanF_maxGFP_noB = ~BoundariesBW.*Merged_replacedF_maxGFP;
    [FTL_tracked_meanF_maxGFP_boundaries] = Merge8bRGB(FTL_tracked_meanF_maxGFP_noB, Boundaries_RGB,show);
    WriteRGB(FTL_tracked_meanF_maxGFP_boundaries, PathToSave, '_maxF_maxGFP+bound.tiff','none')
    clear Boundaries_RGB
    clear FTL_tracked_meanF_maxGFP_noB

    Merged_meanF_maxGFP_RGB(:,:,1,:) = Merged_replacedF_maxGFP(:,:,:);
    Merged_meanF_maxGFP_RGB(:,:,2,:) = Merged_replacedF_maxGFP(:,:,:);
    Merged_meanF_maxGFP_RGB(:,:,3,:) = Merged_replacedF_maxGFP(:,:,:);

    Factor = 2; % 1 in macbook, 2 in pro
    printLabels_new(Merged_meanF_maxGFP_RGB,Stats_Trx,Factor,'off', PathToSave, '_segmented_tracked_info.tiff','packbits')
end

%% remove from here???
%minNumb = 10; Smooth = 3; minOn = 5; minPDis = 10; SplitShortLong = 60; 
%[PropAll] = AnalyzeTraces(0,0,Spots,Stats_GFP,TimeRes,Frames,Slices,Bits,XYRes, ZRes,Width, Height,minNumb,Smooth, minOn,minPDis,SplitShortLong,SplitEarly,nc14,Delay,PathToSave,Nickname,Rep);

disp('finished!')
%
clear variables
close all

end

