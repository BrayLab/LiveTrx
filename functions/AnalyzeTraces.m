%% 
function[Properties] = AnalyzeTraces(varargin)
    manual = varargin{1};
    regions = varargin{2};
    Spots = varargin{3};
    WhichAT = 'MCP';
    try WhichAT = varargin{5}; end
    %% if manual do manual = 1
    %info = readtable('~/Google Drive/MATLAB_R_scripts/metadata MS2 3D.txt','Delimiter','\t')
    if manual
        Parameters = varargin{4};
        %Parameters = info(10,:); % ############
        %minNumb = 5 % ############
        %Spots = 0; % ############
        %regions = 1; % ############
        %x = varargin{5};
        %%
        %Parameters = info(x,:);
        Table2Vars(Parameters);
        Flip = str2double(strsplit(Flip,','));
        PathToSave = [Path,File,Name,File];
        try
            load([PathToSave, '_Stats_tracked.mat']);
        catch
            load([PathToSave, '_Stats.mat']);
            Stats_tracked = Stats_GFP;
            clear Stats_GFP
        end
        try
            SettingsOverwrite = varargin{6};
            cellfun(@(x) evalin('caller',x), SettingsOverwrite);
            %eval(SettingsOverwrite)
            fileID = fopen([PathToSave,'_Settings_AnalyzeTrx_',WhichAT,'.txt'],'a')
            fprintf(fileID, '\nupdated on \t%s . used settings:\n',date)
            fprintf(fileID,'%s \n',SettingsOverwrite{:})
            fclose(fileID);
            SplitShortLong = 60; % ignore
            disp('loaded settings')
        catch
            minNumb = 10; %minimum number of tracked frames for nuclei to be kept
            Smooth = 3; 
            minOn = 5; %minimum number of tracked nuclei is considered actively transcribing to be kept
            FoldOverBaseline = 1.2; %MedFilt > FoldOverBaseline * Baseline  = ON;
            minPDis = 10; %combination of time and space to fix tracking based on sudden ends or starts of transcription
            %minPDis = 15; %#####################
            SplitShortLong = 60; % ignore
        end
        %SplitEarly = 15;
%         try
%             Im = imread([PathToSave, '_MidmRNAProd_2K_noBG.tiff']);
%         catch
%             Im = Read3d([PathToSave, '_MidmRNAProd_2K_noBG.tiff']); Im = Im(:,:,end);
%         end
        F = max(round(35*60/TimeRes)+nc14-Delay,1);
        if F >= Frames ; F = Frames-1; end
        %Width = size(Im,1); Height = size(Im,2);
            try
                Merged_meanF_maxGFP = Read3d([PathToSave, '_replacedF_maxGFP_MCP.tiff']);
            catch
                try
                    Merged_meanF_maxGFP = Read3d([PathToSave, '_replacedF_maxGFP_mean.tiff']);
                catch
                    Merged_meanF_maxGFP = Read3d([PathToSave, '_maxF_maxGFP.tiff']);
                end
            end
            Im = Merged_meanF_maxGFP(:,:,F)./255;
            Width = size(Im,1); 
            Height = size(Im,2);
        if regions %always define regions w MCP. can be reused to analyze other channels
            %boundariesL = Read3d([PathToSave, '_boundariesL.tiff']);
            % rescue boundariesL from StatsGFP
            cmap = jet(1000000);
            cmap_shuffled = cmap(randperm(size(cmap,1)),:);
            boundariesL = zeros(Height,Width,Frames);
            for f = 1:length(Stats_tracked)
                 disp(num2str(f))
                 SingleFrame = Stats_tracked{f};
                 Tracked3D = zeros(Height,Width,Slices);
                 for n = 1:size(SingleFrame,1)
                     Tracked3D(SingleFrame.PixelIdxList{n}) =  SingleFrame.Label(n);
                 end
                 %toTrack(:,:,:,f) = Tracked3D;
                 TrackedBoundL = MAX_proj_3D(Tracked3D);
                 [BoundariesBW, ~, ~] = BoundariesTracked_3D(Tracked3D,cmap_shuffled,'off');
                 boundariesL(:,:,f) = BoundariesBW .* TrackedBoundL;
            end

        end
    else
        [Stats_tracked,TimeRes,Frames,Slices,Bits,XYRes, ZRes,Width, Height,...
            minNumb,Smooth, minOn,minPDis,SplitShortLong,SplitEarly,...
            nc14,Delay,PathToSave,Nickname,Rep] = varargin{7:end}; 
    end
    %%
    % Stats_GFP = Stats_toFix ;
    %%
    try
        Data = Data.Data; % to force error
        if strcmp(WhichAT,'MCP')
            load([PathToSave,'_Data.mat']);
        elseif strcmp(WhichAT,'PCP')
            load([PathToSave,'_PCP_Data.mat']);
        elseif strcmp(WhichAT,'mean')
            load([PathToSave,'_mean_Data.mat']);
        end
        Struct2Vars(Data);
    catch
        %%
        if strcmp(WhichAT,'MCP')
            try
            	load([PathToSave, '_StatsMCP.mat']);
            	[AllF] = MergeAll(Stats_Trx, TimeRes);
            catch
            	[AllF] = MergeAll(Stats_tracked, TimeRes);
            end
        elseif strcmp(WhichAT,'PCP')
                load([PathToSave, '_StatsPCP.mat']);
            	[AllF] = MergeAll(Stats_Trx, TimeRes);
        elseif strcmp(WhichAT,'mean')
                load([PathToSave, '_Statsmean.mat']);
                [AllF] = MergeAll(Stats_Trx, TimeRes);
                PathToSave = [PathToSave,'_mean'];
        end
        AllF = splitvars(AllF);
        %
        Labels = unique(AllF.Label); LabelsOld = Labels;
        [MaxF] = Reshape(AllF,Frames,Labels,'Max','Label');
        [MeanF] = Reshape(AllF,Frames,Labels,'MeanIntensity','Label');
        % MANUAL plot mean nuclear levels
%         Fig1 = figure('PaperSize',[30 50],'PaperUnits','inches','resize','on', 'visible','on');
%         Fig1.Renderer='Painters';
%         [MeanF] = Reshape(AllF,Frames,Labels,'MeanIntensity','Label');
%         [Area] = Reshape(AllF,Frames,Labels,'Area','Label').*XYRes.*XYRes.*ZRes;
%         [CentY] = Reshape(AllF,Frames,Labels,'Centroid_2','Label').*XYRes;
%         TimeScale = ([1: Frames] - nc14+Delay).*TimeRes./60;
%         subplot(311); plot(TimeScale,MeanF,'Color',[0.7,0.7,0.7,0.5]);hold on; plot(TimeScale,nanmean(MeanF,2))
%         xlabel('time into nc14 (min)'); ylabel('Mean Fluorescence (AU)')
%         subplot(312); plot(TimeScale,Area,'Color',[0.7,0.7,0.7,0.5]);hold on; plot(TimeScale,nanmean(Area,2))
%         xlabel('time into nc14 (min)'); ylabel('Volume (um^3)')
%         subplot(313); plot(TimeScale,CentY,'Color',[0.7,0.7,0.7,0.5]);hold on; plot(TimeScale,nanmean(CentY,2))
%         xlabel('time into nc14 (min)'); ylabel('Y Position (um)')
%         print(Fig1,[Path,File,Name,File,'_Ch1.pdf'],'-fillpage', '-dpdf');
%         minNumb = 5 % ############

        % end comment
        
        noNaNnum = sum(~isnan(MaxF),1);
        ToDelete = noNaNnum <= minNumb;
        % delete tracks less than minNumb
        Labels(noNaNnum <= minNumb) = [];
        MaxF(:,noNaNnum <= minNumb) = [];
        MeanF(:,noNaNnum <= minNumb) = [];

        if strcmp(WhichAT,'mean')
            MeasuredF = MeanF; %#########################
        else
            MeasuredF = MaxF;
        end
        
        [Baseline, baseline,MaxFBG, BG, BGLabels] = BaselinefromBG(MeasuredF,FoldOverBaseline,Labels, Smooth,minOn);
        
%         if strcmp(WhichAT,'mean')
%             BG = zeros(size(BG)); % ###########
%             Baseline = Baseline .* 0; % ############
%         end
        % remove BG Labels
        Labels = Labels(~BG);
        %parpool(4)

        %spmd
        tic
        MaxF = CleanMatrixValues(AllF,Frames,LabelsOld,'Max',ToDelete,BG);
        MeanF = CleanMatrixValues(AllF,Frames,LabelsOld,'MeanIntensity',ToDelete,BG);
        Area = CleanMatrixValues(AllF,Frames,LabelsOld,'Area',ToDelete,BG);
        CentX = CleanMatrixValues(AllF,Frames,LabelsOld,'Centroid_1',ToDelete,BG);
        CentY = CleanMatrixValues(AllF,Frames,LabelsOld,'Centroid_2',ToDelete,BG);
        CentZ = CleanMatrixValues(AllF,Frames,LabelsOld,'Centroid_3',ToDelete,BG);
        SpotX = CleanMatrixValues(AllF,Frames,LabelsOld,'SpotPositionX',ToDelete,BG);
        SpotY = CleanMatrixValues(AllF,Frames,LabelsOld,'SpotPositionY',ToDelete,BG);
        SpotZ = CleanMatrixValues(AllF,Frames,LabelsOld,'SpotPositionZ',ToDelete,BG);
        if strcmp(WhichAT,'mean')
            MeasuredF = MeanF; %#########################
        else
            MeasuredF = MaxF;
        end
        % ############
%         Fig1 = figure('PaperSize',[30 50],'PaperUnits','inches','resize','on', 'visible','on');
%         Fig1.Renderer='Painters';
%         subplot(311); plot(TimeScale,MeanF,'Color',[0.7,0.7,0.7,0.5]);hold on; plot(TimeScale,nanmean(MeanF,2))
%         xlabel('time into nc14 (min)'); ylabel('Mean Fluorescence (AU)')
%         subplot(312); plot(TimeScale,Area,'Color',[0.7,0.7,0.7,0.5]);hold on; plot(TimeScale,nanmean(Area,2))
%         xlabel('time into nc14 (min)'); ylabel('Volume (um^3)')
%         subplot(313); plot(TimeScale,CentY,'Color',[0.7,0.7,0.7,0.5]);hold on; plot(TimeScale,nanmean(CentY,2))
%         xlabel('time into nc14 (min)'); ylabel('Y Position (um)')
%         print(Fig1,[Path,File,Name,File,'_Ch1_noNaN.pdf'],'-fillpage', '-dpdf');
        % ############

        if Spots > 0
             [MaxFGauss] = CleanMatrixValues(AllF,Frames,LabelsOld,'MaxGauss',ToDelete,BG);
             [SpotXGauss] = CleanMatrixValues(AllF,Frames,LabelsOld,'PositionGaussX',ToDelete,BG);
             [SpotYGauss] = CleanMatrixValues(AllF,Frames,LabelsOld,'PositionGaussY',ToDelete,BG);
             [SpotZGauss] = CleanMatrixValues(AllF,Frames,LabelsOld,'PositionGaussZ',ToDelete,BG);
             if Spots == 2
                [MaxFGauss2] = CleanMatrixValues(AllF,Frames,LabelsOld,'MaxGauss2',ToDelete,BG);
                [SpotXGauss2] = CleanMatrixValues(AllF,Frames,LabelsOld,'Position2GaussX',ToDelete,BG);
                [SpotYGauss2] = CleanMatrixValues(AllF,Frames,LabelsOld,'Position2GaussY',ToDelete,BG);
                [SpotZGauss2] = CleanMatrixValues(AllF,Frames,LabelsOld,'Position2GaussZ',ToDelete,BG);
             end
        end
        disp('matrices done')
        toc
        %end
 
        %% fix tracking
        %minPDis = 10;
        % modified so that only considers starts and ends during nc14
        [LabelsE,LabelsS] = FixTrackingEnds(MeasuredF,Labels,TimeRes,CentX,CentY,CentZ,minPDis,nc14);
        % replace old and new labels
        %spmd
        [~, MaxF] = ReplaceLabels(Labels,LabelsE,LabelsS,MaxF);
        [~, MeanF] = ReplaceLabels(Labels,LabelsE,LabelsS,MeanF);
        [~, Area] = ReplaceLabels(Labels,LabelsE,LabelsS,Area);
        [~, CentX] = ReplaceLabels(Labels,LabelsE,LabelsS,CentX);
        [~, CentY] = ReplaceLabels(Labels,LabelsE,LabelsS,CentY);
        [~, CentZ] = ReplaceLabels(Labels,LabelsE,LabelsS,CentZ);
        [~, SpotX] = ReplaceLabels(Labels,LabelsE,LabelsS,SpotX);
        [~, SpotY] = ReplaceLabels(Labels,LabelsE,LabelsS,SpotY);
        [~, MeasuredF] = ReplaceLabels(Labels,LabelsE,LabelsS,MeasuredF);
        if Spots > 0
            [~, MaxFGauss] = ReplaceLabels(Labels,LabelsE,LabelsS,MaxFGauss);
            [~, SpotXGauss] = ReplaceLabels(Labels,LabelsE,LabelsS,SpotXGauss);
            [~, SpotYGauss] = ReplaceLabels(Labels,LabelsE,LabelsS,SpotYGauss);
            [~, SpotZGauss] = ReplaceLabels(Labels,LabelsE,LabelsS,SpotZGauss);
            if Spots == 2
                [~,MaxFGauss2] = ReplaceLabels(Labels,LabelsE,LabelsS,MaxFGauss2);
                [~,SpotXGauss2] = ReplaceLabels(Labels,LabelsE,LabelsS,SpotXGauss2);
                [~,SpotYGauss2] = ReplaceLabels(Labels,LabelsE,LabelsS,SpotYGauss2);
                [~,SpotZGauss2] = ReplaceLabels(Labels,LabelsE,LabelsS,SpotZGauss2);
             end
        end
        %subplot(311);Fig;
        %end
        [Labels, SpotZ] = ReplaceLabels(Labels,LabelsE,LabelsS,SpotZ);
        % defining before and after split

        [MedFilt, OnOff,Properties, Fig] = DefineExpAll(MeasuredF,MaxFBG,CentX,Labels,Baseline,FoldOverBaseline,TimeRes,Smooth,SplitShortLong,SplitEarly,nc14,Delay,minOn);
        [Properties] = ReplaceLabelsProps(Labels,LabelsE,LabelsS,Properties);
%         if strcmp(WhichAT,'mean')
%             Properties.Type = string(repmat('ShortMidline',height(Properties),1)) % ########################
%         end        
        disp('tracking fixed')
        % save MaxF and Properties, done
        % add split exp to file
        %%
        TimeScale = ([1: Frames] - nc14+Delay).*TimeRes./60;
        Data = struct('MaxF', MaxF,'MedFilt', MedFilt,'OnOff',OnOff,...
            'CentX',CentX,'CentY',CentY,'CentZ',CentZ,'SpotX',SpotX,'SpotY',SpotY,'SpotZ',SpotZ,...
            'Baseline',Baseline, 'TimeScale',TimeScale,...
            'Properties',Properties,'MaxFBG',MaxFBG,'BGLabels',BGLabels,'MeanF',MeanF,'Area',Area);
        if Spots > 0
            Data.MaxFGauss = MaxFGauss;
            Data.SpotXGauss = SpotXGauss;
            Data.SpotYGauss = SpotYGauss;
            Data.SpotZGauss = SpotZGauss;
            if Spots == 2
                Data.MaxFGauss2 = MaxFGauss2;
                Data.SpotXGauss2 = SpotXGauss2;
                Data.SpotYGauss2 = SpotYGauss2;
                Data.SpotZGauss2 = SpotZGauss2;
            end
        end
        save([PathToSave,'_Data.mat'],'Data');
        

    end
    %% PLOTS
    %plot selected
    %set(0,'defaultAxesFontSize',14)
    %set(0,'defaultAxesFontName', 'Myriad Pro')
    
    Selected = [Properties.Type=='EarlyOnly'|Properties.Type=='ShortMidline'|Properties.Type=='LongMidline'];
    if sum(Selected) == 0
        Selected = logical(Selected +1);
    end
    % count bursts and frequency
    [BurstNum,BurstLength,BurstPeriod,BurstPeak,BurstMax,OffTimeAll,BurstSize] = CountBursts(MeasuredF,OnOff, Selected,minOn,1,TimeRes);
    TimeScale = ([1: Frames] - nc14+Delay).*TimeRes./60;
    Limits = [min((1 - nc14+Delay).*TimeRes./60,1), 90];
    
    PlotSelected(Selected, TimeScale,Baseline,MeasuredF,MedFilt,OnOff,Properties,BurstPeak,TimeRes,Bits,nc14, Delay,Limits,[PathToSave,'_selected.ps']);
    if Spots > 0
        PlotSelected(Selected, TimeScale,Baseline,MaxFGauss,MeasuredF,OnOff,Properties,BurstPeak,TimeRes,Bits+2,nc14, Delay,Limits,[PathToSave,'_selected_spots.ps']);
        if Spots == 2
            PlotSelected(Selected, TimeScale,Baseline,MaxFGauss,MaxFGauss2,OnOff,Properties,BurstPeak,TimeRes,Bits+2,nc14, Delay,Limits,[PathToSave,'_selected_2spots.ps']);
        end
    end
    
    if ~strcmp(WhichAT,'mean')
    try
    % plot spot position relative to centroid for all midline
    Selected = [Properties.Type=='ShortMidline'|Properties.Type=='LongMidline'];
    SpotXOn = SpotX .* OnOff; SpotXOn(SpotXOn==0) = NaN;
    SpotYOn = SpotY .* OnOff; SpotYOn(SpotYOn==0) = NaN;
    if Spots == 0
        PlotSpotMov([PathToSave,'_SpotsMov.ps'], Spots, Selected,Properties, CentX, CentY, SpotXOn, SpotYOn)
    end
    if Spots == 1
        PlotSpotMov([PathToSave,'_SpotsMov.ps'], Spots, Selected,Properties, CentX, CentY, SpotXOn, SpotYOn,SpotXGauss,SpotYGauss)
    end
    if Spots == 2
        PlotSpotMov([PathToSave,'_SpotsMov.ps'], Spots, Selected,Properties, CentX, CentY, SpotXOn, SpotYOn,SpotXGauss,SpotYGauss,SpotXGauss2,SpotYGauss2)
    end
    % so that it prints properly, Selected has to be the same used to count
    % bursts
    end
    end
    %%
    MaxTracked = zeros(Width,Height,Frames);
    for f = 1:length(Stats_tracked)
         SingleFrame = Stats_tracked{f};
         Tracked3D = zeros(Width,Height,Slices);
         for n = 1:size(SingleFrame,1)
             Tracked3D(SingleFrame.PixelIdxList{n}) =  SingleFrame.Label(n);
         end
         MaxTracked(:,:,f) = max(Tracked3D,[],3);
    end
    %%
    Selected = [Properties.Type ~= 'BG'];
    %MaxFOn = MaxF.*OnOff; MaxFOn(MaxFOn==0) = NaN; 
    Norm = (MeasuredF-Baseline').*Baseline(1)./Baseline'.*2^(12-Bits);
    NormSelected = Norm(:,Selected);
    NormBGB = (MaxFBG-Baseline').*Baseline(1)./Baseline'.*2^(12-Bits);
    PropertiesSelected = Properties(Selected,:);
    ElT = 5000/2000*60;
    SplitEarlyF = max(floor([SplitEarly*60./TimeRes+nc14-Delay,1]));
    TotalmRNAMid = nansum(NormSelected(SplitEarlyF+1:end,:),1) * TimeRes ./ 2 ./ ElT;
    TotalmRNAEarly = nansum(NormSelected(1:SplitEarlyF,:),1) * TimeRes ./ 2 ./ ElT;
    TotalmRNABGMid = nansum(NormBGB(SplitEarlyF+1:end,:),1) * TimeRes ./ 2 ./ ElT;
    TotalmRNABGEarly = nansum(NormBGB(1:SplitEarlyF,:),1) * TimeRes ./ 2 ./ ElT;
    %plot(PropertiesSelected.AP_position,TotalmRNAMid,'o')
    %plot(PropertiesSelected.TimeOn,TotalmRNAMid,'o')
    Properties.TotalmRNA = TotalmRNAMid';
      Data.Properties = Properties;
              save([PathToSave,'_Data.mat'],'Data');
    
    SelectedLabels = PropertiesSelected.Label;
    %MaxTrackedMid = MaxTracked(:,:,SplitEarlyF+1:end);
    %MaxTrackedEarly = MaxTracked(:,:,1:SplitEarlyF);
    MaxTrackedMid = MaxTracked(:,:,F); %for single frame early
    MaxTrackedEarly = MaxTracked(:,:,SplitEarlyF); % for single frame mid
    TotalProdMid = zeros(size(MaxTrackedMid));
    TotalProdEarly = zeros(size(MaxTrackedEarly));
    
    for n = 1:length(SelectedLabels)
        %disp([num2str(n),'/',num2str(length(SelectedLabels)),' ON'])
        NewLabels = PropertiesSelected.NewLabel(n);
        TotalProdMid(ismember(MaxTrackedMid,NewLabels{:})) =  TotalmRNAMid(n);
        TotalProdEarly(ismember(MaxTrackedEarly,NewLabels{:})) =  TotalmRNAEarly(n);
    end
    for n = 1:length(BGLabels)
        %disp([num2str(n),'/',num2str(length(BGLabels)),' BG'])
        TotalProdMid(ismember(MaxTrackedMid,BGLabels(n))) =  1; %noBG
        TotalProdEarly(ismember(MaxTrackedEarly,BGLabels(n))) =  1; %noBG
        %TotalProdMid(ismember(MaxTrackedMid,BGLabels(n))) =  TotalmRNABGMid(n);
        %TotalProdEarly(ismember(MaxTrackedEarly,BGLabels(n))) =  TotalmRNABGEarly(n);
    end
    %
    maxTotRNA = 10000; % 3000 for lateral. 10000 for ventral
    %TotalProd = cat(3,TotalProdEarly,TotalProdMid);
    TotalProd = TotalProdMid; % for single frame
    TotalProd (TotalProd < 0) = 1;
    TotalProd (TotalProd > maxTotRNA) = maxTotRNA;
    TotalProdRBG = zeros(size(TotalProd,1), size(TotalProd,2),3,size(TotalProd,3));
    CMAP = colormap(parula(maxTotRNA+1));
    for f = 1:size(TotalProd,3)
        TotalProdRBG(:,:,:,f) = label2rgb(round(TotalProd(:,:,f)),CMAP,'k');
    end
    WriteRGB(double(TotalProdRBG), PathToSave, '_MidmRNAProd_10K_noBG.tiff','none')
    
    maxTotRNA = 2000; % 3000 for lateral. 10000 for ventral
    %TotalProd = cat(3,TotalProdEarly,TotalProdMid);
    TotalProd = TotalProdMid;
    TotalProd (TotalProd < 0) = 1;
    TotalProd (TotalProd > maxTotRNA) = maxTotRNA;
    TotalProdRBG = zeros(size(TotalProd,1), size(TotalProd,2),3,size(TotalProd,3));
    CMAP = colormap(parula(maxTotRNA+1));   
    for f = 1:size(TotalProd,3)
        TotalProdRBG(:,:,:,f) = label2rgb(round(TotalProd(:,:,f)),CMAP,'k');
    end
    WriteRGB(double(TotalProdRBG), PathToSave, '_MidmRNAProd_2K_noBG.tiff','none')
  

     %%

    %% amplitude, frequency and length
    %
    % count bursts and frequency
    Selected = [Properties.Type=='ShortMidline'|Properties.Type=='LongMidline'];
    if sum(Selected) == 0
        Selected = logical(Selected +1);
    end
    
    SplitEarlyF = max([SplitEarly*60./TimeRes+nc14-Delay,1]);
    [BurstNum,BurstLength,BurstPeriod,BurstPeak,BurstMax,OffTimeAll, BurstSize] = CountBursts(MeasuredF,OnOff, Selected,minOn,SplitEarlyF,TimeRes);
    PlotBurstingProp(BurstNum,BurstLength,BurstPeriod,TimeRes,[PathToSave,'_selected_plots.ps'])
    
     Selected = [(Properties.Type=='ShortMidline'|Properties.Type=='LongMidline') & Properties.TimeOn > 5];
     if sum(Selected) == 0
        Selected = [Properties.Type~='EarlyOnly'];
    end
%      PropertiesSelected = Properties (Selected,:);
%      subplot(131);plot([rand(1,length(PropertiesSelected.Onset))],PropertiesSelected.Onset,'.r'); hold on
%      scatter([rand(1,length(PropertiesSelected.End))],PropertiesSelected.End,'.b'); ylim([0,90])
%      subplot(132);boxplot([PropertiesSelected.Onset,PropertiesSelected.End],'Jitter',0.5,'Symbol','.','Whisker',1); ylim([0,90])
%      subplot(133);plot(PropertiesSelected.AP_position,PropertiesSelected.Onset,'.r'); hold on
%      plot(PropertiesSelected.AP_position,PropertiesSelected.End,'.b');  ylim([0,90])

    % 
    % ColorArg is 8 value colormap with info for allPoints, allMidPoints,
    % shortPoints, longPoints, allMeans, allMidMeans, shortMeans, longMeans
    % in default it is using all/short/long + grey for all mid, all means
    % red
    cmap_5 = parula(10); %2, 6, 9 for earlyonly, short, long
    ColorArg = [cmap_5(2,:);0.7,0.7,0.7;cmap_5(6,:);cmap_5(9,:);...
        1,0,0;1,0,0;1,0,0;1,0,0];    
    Fig = figure('PaperSize',[50 50],'PaperUnits','inches','resize','on', 'visible','off');
    Fig.Renderer='Painters';
    Fig = PlotMeans(MeasuredF,TimeScale,Selected,Properties,Bits,ColorArg,Fig,0,[Nickname,' ',num2str(Rep)],[-50,90], [min(0,min(MeasuredF(:))), max(2^Bits-1,max(MeasuredF(:)))]);
    print(Fig,[PathToSave,'_selected_plots.ps'],'-fillpage', '-dpsc','-append');
    %
    PlotsPosition(OnOff,TimeScale,Selected,Properties,CentX,CentY,CentZ,SpotZ,XYRes,ZRes,Width, Height,[PathToSave,'_selected_plots.ps'])
    %%
    Im = imread([PathToSave, '_MidmRNAProd_10K_noBG.tiff']);

    if regions
         mkdir([Path,File,Name,'regions/'])
         PathToSave = [Path,File,Name,'regions/',File]; 
         close all
        %%
        try
            Regions = imread([PathToSave,'_regions.tiff']);
            if regions == 2
                error('overwritting regions.tiff')
            end
            %Old = imread([PathToSave,'_regions_old.tiff']);
            %imwrite(uint8(Regions), [PathToSave,'_regions_old.tiff'],'Compression','none')
            cmap_3 = [255,120,120;096,085,176;95,181,241;80,159,115]./255;
            mkdir([Path,File,Name,'regions/'])
            PathToSave = [Path,File,Name,'regions/',File]; 
            close all
        catch
            figure;plot(CentX,-CentY); pause(1);
%             SubIm = Im(:,:,1);
            %[CentY(:)]
%             Indices = [sub2ind([size(SubIm)],[round(CentY(:))],[round(CentX(:))])]';
%             Indices(isnan(Indices)) = []
%             SubIm(Indices) = 1;
%             Im(:,:,1) = SubIm;
            figure;
            RME = roipoly(Im);
            RMSE = roipoly(Im.*repmat(uint8(~RME),[1,1,3]));RMSE = RMSE.*2;
            RNE = roipoly(Im.*repmat(uint8(~(RME|RMSE)),[1,1,3])); RNE = RNE.*3;
            RDE = roipoly(Im.*repmat(uint8(~(RME|RMSE|RNE)),[1,1,3])); RDE = RDE.*4;
            %Regions = RME + RMSE + RNE + RDE;
            Regions = cat(3,RME,RMSE,RNE,RDE);
            Regions = max(Regions,[],3);
            ToDraw = max(cat(3,bwmorph(RME,'remove'),bwmorph(RMSE,'remove').*2,bwmorph(RNE,'remove').*3,bwmorph(RDE,'remove').*4),[],3);
            cmap_3 = [255,120,120;096,085,176;95,181,241;80,159,115]./255;
            ToDraw = label2rgb(ToDraw, cmap_3,'k');
            imagesc(uint8(Im*255)); hold on
            Merged = Merge8bRGB(uint8(Im*255), ToDraw,'off');
            imagesc(Merged); colormap(cmap_3)
            imwrite(uint8(Regions), [PathToSave,'_regions.tiff'],'Compression','none')
            imwrite(Merged, [PathToSave,'_regionsMerged.tiff'],'Compression','none')
        end
        if strcmp(WhichAT,'mean')
            PathToSave = [Path,File,Name,'regions/',File,'_mean'];
        else
            PathToSave = [Path,File,Name,'regions/',File]; 
        end
        %imagesc(Im,[0,1])    
        %% find centroids of each region in selected frame
        Selected = find([Properties.Type=='ShortMidline'|Properties.Type=='LongMidline']);
        if sum(Selected) == 0 %only when movie stars later then get LateTracks
            Selected = [Properties.Type~='EarlyOnly'];
         end
        %LabelsSelected = Labels(Selected);
        PosX = floor(nanmean(CentX(max(1,F-5): min(F+5,Frames-1),Selected),1));
        PosY = floor(nanmean(CentY(max(1,F-5): min(F+5,Frames-1),Selected),1));
        %Indices = sub2ind(size(Regions),PosY,PosX);
        % added to prevent too many NaNs when cells are not tracked in the
        % F-5:F+5 window
        PosX(isnan(PosX)) = floor(nanmean(CentX(max(1,nc14): min(Frames-1,nc14-Delay+40*60./TimeRes),Selected(isnan(PosX))),1));
        PosY(isnan(PosY)) = floor(nanmean(CentY(max(1,nc14): min(Frames-1,nc14-Delay+40*60./TimeRes),Selected(isnan(PosY))),1));
        Indices = sub2ind(size(Regions),PosY,PosX);
        Selected = Selected(~isnan(Indices));
        %LabelsSelected = LabelsSelected(~isnan(Indices));
        PosX = PosX(~isnan(Indices));
        PosY = PosY(~isnan(Indices));
        Indices = Indices(~isnan(Indices));
        RegionsInd = Regions(Indices);
        ME = Selected(RegionsInd==1);
        MSE = Selected(RegionsInd==2);
        NE = Selected(RegionsInd==3);
        DE = Selected(RegionsInd==4);
        Properties.Region = string(repmat('NaN',length(Properties.Type),1));
        Properties.Region(ME) = string(repmat('ME',length(ME),1));
        Properties.Region(MSE) = string(repmat('MSE',length(MSE),1));
        Properties.Region(NE) = string(repmat('NE',length(NE),1));
        Properties.Region(DE) = string(repmat('DE',length(DE),1));
        %% fit line to MSE cells (if none it is horizontal height = 0 bar) and calculate distance from each to MSE
        %Properties.DMSE = [];
        XcoordMSE = PosX(RegionsInd==2);
        YcoordMSE = PosY(RegionsInd==2);
        if isempty(XcoordMSE)
            XcoordMSE = 1:1:Width;
            YcoordMSE = zeros(1,Width)+Height;
        end
        p = polyfit(XcoordMSE,YcoordMSE,2);
        a = p(1); b = p(2); c = p(3);
        MinimumDistance = arrayfun(@(X,Y) max(roots([2*c^2, 3*a*b, b^2+2*a*c-2*c*Y+1, a*b-b*Y-X])), PosX,PosY,'UniformOutput', false);
        MinimumDistance = cell2mat(MinimumDistance);
        Properties.DMSE(Selected) = MinimumDistance;
        Properties.DMSE(Properties.DMSE==0) = NaN;
        %%
        % ColorArg is 8 value colormap with info for allPoints, allMidPoints,
        % shortPoints, longPoints, allMeans, allMidMeans, shortMeans, longMeans
        % in default it is using all/short/long + grey for all mid, all means
        % red
        % for regions 1-4 values are replaced with RGB for each region,
        % means still red
        TimeScale = ([1:Frames] - nc14+Delay).*TimeRes./60;
        Limits = [0, 90];
        try
        Selected = ME; 
        VarName = varname(ME);
        ColorArg = [cmap_3(1,:);cmap_3(1,:);cmap_3(1,:);cmap_3(1,:);...
            1,0,0;1,0,0;1,0,0;1,0,0];    
        [BurstNum,BurstLength,BurstPeriod,BurstPeak,BurstMax,OffTimeAll, BurstSize] = CountBursts(MeasuredF,OnOff, Selected,minOn,SplitEarlyF,TimeRes);
        PlotSelected(Selected, TimeScale,Baseline,MeasuredF,MedFilt,OnOff,Properties,BurstPeak,TimeRes,Bits,nc14, Delay,Limits,[PathToSave,'_selected_',VarName,'.ps']);
        PlotBurstingProp(BurstNum,BurstLength,BurstPeriod,TimeRes,[PathToSave,'_selected_plots_',VarName,'.ps'])
        Fig = figure('PaperSize',[50 50],'PaperUnits','inches','resize','on', 'visible','off');
        Fig = PlotMeans(MeasuredF,TimeScale,Selected,Properties,Bits,ColorArg,Fig,0,[Nickname,' ',num2str(Rep),' ',VarName],[-50,90], [min(0,min(MeasuredF(:))), max(2^Bits-1,max(MeasuredF(:)))]);
        print(Fig,[PathToSave,'_selected_plots_',VarName,'.ps'],'-fillpage', '-dpsc','-append');
        PlotsPosition(OnOff,TimeScale,Selected,Properties,CentX,CentY,CentZ,SpotZ,XYRes,ZRes,Width, Height,[PathToSave,'_selected_plots_',VarName,'.ps'])
        end
        try
        Selected = MSE; 
        VarName = varname(MSE);
        ColorArg = [cmap_3(2,:);cmap_3(2,:);cmap_3(2,:);cmap_3(2,:);...
            1,0,0;1,0,0;1,0,0;1,0,0];
        [BurstNum,BurstLength,BurstPeriod,BurstPeak,BurstMax,OffTimeAll,BurstSize] = CountBursts(MeasuredF,OnOff, Selected,minOn,SplitEarlyF,TimeRes);
        PlotSelected(Selected, TimeScale,Baseline,MeasuredF,MedFilt,OnOff,Properties,BurstPeak,TimeRes,Bits,nc14, Delay,Limits,[PathToSave,'_selected_',VarName,'.ps']);
        PlotBurstingProp(BurstNum,BurstLength,BurstPeriod,TimeRes,[PathToSave,'_selected_plots_',VarName,'.ps'])
        Fig = figure('PaperSize',[50 50],'PaperUnits','inches','resize','on', 'visible','off');
        Fig = PlotMeans(MeasuredF,TimeScale,Selected,Properties,Bits,ColorArg,Fig,0,[Nickname,' ',num2str(Rep),' ',VarName],[-50,90], [min(0,min(MeasuredF(:))), max(2^Bits-1,max(MeasuredF(:)))]);
        print(Fig,[PathToSave,'_selected_plots_',VarName,'.ps'],'-fillpage', '-dpsc','-append');
        PlotsPosition(OnOff,TimeScale,Selected,Properties,CentX,CentY,CentZ,SpotZ,XYRes,ZRes,Width, Height,[PathToSave,'_selected_plots_',VarName,'.ps'])
        end
        try
        Selected = NE; 
        VarName = varname(NE);
        ColorArg = [cmap_3(3,:);cmap_3(3,:);cmap_3(3,:);cmap_3(3,:);...
            1,0,0;1,0,0;1,0,0;1,0,0];
        [BurstNum,BurstLength,BurstPeriod,BurstPeak,BurstMax,OffTimeAll,BurstSize] = CountBursts(MeasuredF,OnOff, Selected,minOn,SplitEarlyF,TimeRes);
        PlotSelected(Selected, TimeScale,Baseline,MeasuredF,MedFilt,OnOff,Properties,BurstPeak,TimeRes,Bits,nc14, Delay,Limits,[PathToSave,'_selected_',VarName,'.ps']);
        PlotBurstingProp(BurstNum,BurstLength,BurstPeriod,TimeRes,[PathToSave,'_selected_plots_',VarName,'.ps'])
        Fig = figure('PaperSize',[50 50],'PaperUnits','inches','resize','on', 'visible','off');
        Fig = PlotMeans(MeasuredF,TimeScale,Selected,Properties,Bits,ColorArg,Fig,0,[Nickname,' ',num2str(Rep),' ',VarName],[-50,90], [min(0,min(MeasuredF(:))), max(2^Bits-1,max(MeasuredF(:)))]);
        print(Fig,[PathToSave,'_selected_plots_',VarName,'.ps'],'-fillpage', '-dpsc','-append');
        PlotsPosition(OnOff,TimeScale,Selected,Properties,CentX,CentY,CentZ,SpotZ,XYRes,ZRes,Width, Height,[PathToSave,'_selected_plots_',VarName,'.ps'])
        end
        try
        Selected = DE; 
        VarName = varname(DE);
        ColorArg = [cmap_3(4,:);cmap_3(4,:);cmap_3(4,:);cmap_3(4,:);...
            1,0,0;1,0,0;1,0,0;1,0,0];
        [BurstNum,BurstLength,BurstPeriod,BurstPeak,BurstMax,OffTimeAll,BurstSize] = CountBursts(MeasuredF,OnOff, Selected,minOn,SplitEarlyF,TimeRes);
        PlotSelected(Selected, TimeScale,Baseline,MeasuredF,MedFilt,OnOff,Properties,BurstPeak,TimeRes,Bits,nc14, Delay,Limits,[PathToSave,'_selected_',VarName,'.ps']);
        PlotBurstingProp(BurstNum,BurstLength,BurstPeriod,TimeRes,[PathToSave,'_selected_plots_',VarName,'.ps'])
        Fig = figure('PaperSize',[50 50],'PaperUnits','inches','resize','on', 'visible','off');
        Fig = PlotMeans(MeasuredF,TimeScale,Selected,Properties,Bits,ColorArg,Fig,0,[Nickname,' ',num2str(Rep),' ',VarName],[-50,90], [min(0,min(MeasuredF(:))), max(2^Bits-1,max(MeasuredF(:)))]);
        print(Fig,[PathToSave,'_selected_plots_',VarName,'.ps'],'-fillpage', '-dpsc','-append');
        PlotsPosition(OnOff,TimeScale,Selected,Properties,CentX,CentY,CentZ,SpotZ,XYRes,ZRes,Width, Height,[PathToSave,'_selected_plots_',VarName,'.ps'])
        end
        close all
        %%
        MELabels = Properties.NewLabel(ME); MELabels = [MELabels{:}];
        MSELabels = Properties.NewLabel(MSE); MSELabels = [MSELabels{:}];
        NELabels = Properties.NewLabel(NE); NELabels = [NELabels{:}];
        DELabels = Properties.NewLabel(DE); DELabels = [DELabels{:}];

        boundariesL_ME = ismember(boundariesL,MELabels');
        boundariesL_MSE = ismember(boundariesL,MSELabels')*2;
        boundariesL_NE = ismember(boundariesL,NELabels')*3;
        boundariesL_DE = ismember(boundariesL,DELabels')*4;
        boundariesL_selected = boundariesL_ME + boundariesL_MSE + boundariesL_NE + boundariesL_DE;

        boundariesL_selected_RGB = zeros(size(boundariesL_selected,1),size(boundariesL_selected,2),3,size(boundariesL_selected,4));
        for f = [1:size(boundariesL_selected,3)]
            boundariesL_selected_RGB(:,:,:,f) = label2rgb(boundariesL_selected(:,:,f), cmap_3, 'k', 'noshuffle');
        end  
        FTL_tracked_meanF_maxGFP_noB_selected = ~boundariesL_selected.*Merged_meanF_maxGFP;
        [FTL_tracked_meanF_maxGFP_boundaries_selected] = Merge8bRGB(FTL_tracked_meanF_maxGFP_noB_selected, boundariesL_selected_RGB,'off');
        [FTL_tracked_meanF_maxGFP_boundaries_selected] = TimeStamp(uint8(FTL_tracked_meanF_maxGFP_boundaries_selected),TimeRes,nc14,Delay);
        WriteRGB(double(FTL_tracked_meanF_maxGFP_boundaries_selected), PathToSave, '_maxF+SelAll_regions.tiff','none')
        %
        % Factor = 2; % 1 in macbook, 2 in pro
        % printLabels_new(FTL_tracked_meanF_maxGFP_boundaries_selected,Stats_GFP,Factor,'off', PathToSave, '_segmented_tracked_info_regions.tiff','packbits')
       
    %   
        Data.Properties = Properties;
        if strcmp(WhichAT,'mean')
            PathToSave = [Path,File,Name,File,'_mean'];
        else
            PathToSave = [Path,File,Name,File]; 
        end
        save([PathToSave,'_Data.mat'],'Data');
    end
     disp('finished!')
close all
end

