function[] = CompareBursting(info,Nicknames,ExpLabels,Selections,Palette,PathMain,varargin)

% Nicknames = Input{e}.Nicknames;
% ExpLabels = Input{e}.ExpLabels;
% Selections = Input{e}.Selections;
% Palette = Input{e}.Palette;  

% default settings    

% SplitEarly loaded automatically from Table2Vars, add override?
% and bursts counted from SplitEarly:end

try
    Command = varargin{1};
    eval(Command);
catch
    Log = 1;
    minOn = 5;
    MinContinuous = 10; %min t (min) at least one burst length
    Jitter = 0.6; %Jitter./2 cant be > BarW
    BarW = 0.4;
    FaceAlpha = 0.3;
    DotSize = 5;
    LineWidth = 1.5;
    FontSizeTitle = 10;
    FontSize = 8;
end

MaxRow = length(Nicknames);
BNumAll = cell(1,MaxRow);
BLenAll = cell(1,MaxRow);
MaxLenAll = cell(1,MaxRow);
BPerAll = cell(1,MaxRow);
BAmpAll = cell(1,MaxRow);
OffTimeAll = cell(1,MaxRow);
BSizeAll = cell(1,MaxRow);
OnsetAll = cell(1,MaxRow);
EndAll = cell(1,MaxRow);
TmAll = cell(1,MaxRow);


NumOn = [];
PercOn = [];
PercCont = [];
MaxCol = 1;

for Exp = 1:length(Nicknames)
    % select from REGEXP mode
    %Index = find(cellfun(@(x) ~isempty(x),regexp(info.Nickname, Nicknames{Exp},'match')))
    % select from classic matching  exp
    Index = find(cellfun(@(x) strcmp(x,Nicknames{Exp}),info.Nickname)==1)';
    info.Nickname(Index)
     MaxCol = max(MaxCol,length(Index));
    Selection = Selections{Exp};
for i = [1:length(Index)]
        x = Index(i);
        %
        Parameters = info(x,:);
        Table2Vars(Parameters);
        Flip = str2double(strsplit(Flip,','));
        PathToSave = [Path,File,Name,File]; 
        %load([PathToSave, '_Stats.mat']);
        %minNumb = 10; Smooth = 3; minOn = 5; minPDis = 10; SplitShortLong = 60; SplitEarly = 15;
        Merged_meanF_maxGFP = Read3d([PathToSave, '_maxF_maxGFP.tiff']);
        F = max(round(35*60/TimeRes)+nc14-Delay,1);
        Im = Merged_meanF_maxGFP(:,:,F)./255;
        Width = size(Im,1); Height = size(Im,2);
        clear CentroidF;
        load([PathToSave,'_Data.mat']);
        if length(fields(Data)) > 2
            Struct2Vars(Data);
        else
            Data = Data.Data;
            Struct2Vars(Data);
        end
        
        try
            ImLab = CentroidsF;
        catch
            load([PathToSave, '_Stats.mat']);
            Merged_meanF_maxGFP = Read3d([PathToSave, '_maxF_maxGFP.tiff']);
            Im = Merged_meanF_maxGFP(:,:,F)./255;
            Width = size(Im,1); 
            Height = size(Im,2);
             % rescue boundariesL from StatsGFP
            cmap = jet(1000000);
            cmap_shuffled = cmap(randperm(size(cmap,1)),:);
            boundariesL = zeros(Height,Width,Frames);
            for f = 1:length(Stats_GFP)
                 disp(num2str(f))
                 SingleFrame = Stats_GFP{f};
                 Tracked3D = zeros(Height,Width,Slices);
                 for n = 1:size(SingleFrame,1)
                     Tracked3D(SingleFrame.PixelIdxList{n}) =  SingleFrame.Label(n);
                 end
                 %toTrack(:,:,:,f) = Tracked3D;
                 TrackedBoundL = MAX_proj_3D(Tracked3D);
                 [BoundariesBW, ~, ~] = BoundariesTracked_3D(Tracked3D,cmap_shuffled,'off');
                 boundariesL(:,:,f) = BoundariesBW .* TrackedBoundL;
            end   
            
            ImLab = boundariesL(:,:,max(1,F-5):min(F+5,Frames));
            for f = 1:size(ImLab,3)
                Stats{f} = regionprops('table',ImLab(:,:,f),'Centroid');
                Stats{f}.Label = [1:height(Stats{f})]';
            end
            [AllF] = MergeAll(Stats, TimeRes);
            AllF = splitvars(AllF);
            Labels = unique(AllF.Label); LabelsOld = Labels;
            [CentXF] = Reshape(AllF,Frames,Labels,'Centroid_1','Label');
            [CentYF] = Reshape(AllF,Frames,Labels,'Centroid_2','Label');
            PosX = nanmean(CentXF,1);
            PosY = nanmean(CentYF,1);
             PosX(isnan(PosX)) = [];
             PosY(isnan(PosY)) = [];
%             StatsF = Stats_GFP{F};
%             CentXF = StatsF.Centroid(:,1);
%             CentYF = StatsF.Centroid(:,2);
            Indices = sub2ind([size(Im,1),size(Im,2)],floor(PosY), floor(PosX));
            ImLab = zeros(size(Im,1),size(Im,2));
            ImLab(Indices) = 1;
            Data.CentroidsF = ImLab;
            save([PathToSave,'_Data.mat'],'Data');
        end
        


      %

        if strcmp(Selection,'ME') || strcmp(Selection,'MSE') || strcmp(Selection,'NE') || strcmp(Selection,'DE') == 1
            %mkdir([Path,File,Name,'regions/'])
             PathToSave = [Path,File,Name,'regions/',File]; 
            %
                Regions = imread([PathToSave,'_regions.tiff']);
                cmap_3 = [255,120,120;096,085,176;95,181,241]./255;
                mkdir([Path,File,Name,'regions/'])
                PathToSave = [Path,File,Name,'regions/',File]; 
                Selected = find(Properties.Type~='EarlyOnly');
            %LabelsSelected = Labels(Selected);
            PosX = floor(nanmean(CentX(max(1,F-5): min(F+5,Frames),Selected),1));
            PosY = floor(nanmean(CentY(max(1,F-5): min(F+5,Frames),Selected),1));
            Indices = sub2ind(size(Regions),PosY,PosX);
            Selected = Selected(~isnan(Indices));
            %LabelsSelected = LabelsSelected(~isnan(Indices));
            Indices = Indices(~isnan(Indices));
            RegionsInd = Regions(Indices);
            switch Selection
                case 'ME'
                    Value = 1;
                case 'MSE' 
                    Value = 2;
                case 'NE'
                    Value = 3;
                case 'DE' 
                    Value = 4;
            end
            
            Selected = Selected(RegionsInd==Value);
            %TotalCells = length(unique(ImLab.*(Regions == Value))) - 1;
            TotalCells = ImLab.*(Regions == Value);
            TotalCells = sum(TotalCells(:));
            %figure;imagesc(ImLab)
        end
        
        if strcmp(Selection,'EarlyOnly')
            Selected = find(Properties.Type=='EarlyOnly'); 
            TotalCells = length(unique(ImLab)) - 1;
        end
        
        if strcmp(Selection,'')
            Selected = find(Properties.Type ~= 'EarlyOnly'); 
            TotalCells = length(unique(ImLab)) - 1;
        end
        SplitEarlyF = max([SplitEarly*60./TimeRes+nc14-Delay,1]);
        %VarName = varname(ME);
        NormF = (MaxF-Baseline').*Baseline(1)./Baseline'.*2.^(12-Bits);
        OnOff = CleanOnOff(OnOff,minOn);
        [OnOff] = CleanNaNs(MedFilt,OnOff, minOn*2);
        try
            [BurstNum,BurstLength,BurstPeriod,BurstPeak,BurstMax,OffTime,BurstSize] = CountBursts(NormF,OnOff, Selected,minOn,SplitEarlyF,TimeRes);
        catch
            BurstNum = {}; BurstLength = {}; BurstPeriod = {}; BurstPeak = {};BurstMax = {}; OffTime = {}; BurstSize = {};
        end
        Onsets = Properties.Onset(Selected);
        Ends = Properties.End(Selected);
        TotmRNA = Properties.TotalmRNA(Selected);

        % % cells on
        NumOn(i,Exp) = length(Selected);
        PercOn(i,Exp) = length(Selected) ./ TotalCells;
        Continuous = cellfun(@(x)any(x > MinContinuous.*60./TimeRes),BurstLength);
        PercCont(i,Exp) = sum(Continuous)./length(Continuous);
        NumOn(NumOn == 0) = NaN;
        PercCont(PercOn == 0) = NaN;
        PercOn(PercOn == 0) = NaN;
        PercOn(PercOn > 1) = 1;
        MaxLengthEach = cellfun(@(x) max(x), BurstLength,'UniformOutput',false);
        BNumAll{Exp} = [BNumAll{Exp}; [BurstNum{:}]'];
        BLenAll{Exp} = [BLenAll{Exp}; abs([BurstLength{:}]'.*TimeRes./60)];
        MaxLenAll{Exp} = [MaxLenAll{Exp}; abs([MaxLengthEach{:}]'.*TimeRes./60)];
        BPerAll{Exp} = [BPerAll{Exp}; [BurstPeriod{:}]'.*TimeRes./60];
        BAmpAll{Exp} = [BAmpAll{Exp}; [BurstMax{:}]'];
        OffTimeAll{Exp} = [OffTimeAll{Exp}; [OffTime{:}]'.*TimeRes./60];
        BSizeAll{Exp} = [BSizeAll{Exp}; [BurstSize{:}]'];
        OnsetAll{Exp} = [OnsetAll{Exp}; [Onsets(:)]];
        EndAll{Exp} = [EndAll{Exp}; [Ends(:)]];
        TmAll{Exp} = [TmAll{Exp}; [TotmRNA(:)]];



end
end
%
BNumAll = Cell2Mat(BNumAll);
BLenAll = Cell2Mat(BLenAll);
MaxLenAll = Cell2Mat(MaxLenAll);
BPerAll = Cell2Mat(BPerAll);
BAmpAll = Cell2Mat(BAmpAll);
OffTimeAll = Cell2Mat(OffTimeAll);
BSizeAll = Cell2Mat(BSizeAll);
OnsetAll = Cell2Mat(OnsetAll);
EndAll = Cell2Mat(EndAll);
TmAll = Cell2Mat(TmAll);




%%
mkdir(PathMain)
FileOut = [PathMain,'/comp_bursting_violin',cell2mat(join(unique(Selections),'_')),'_',cell2mat(join(Nicknames,'_vs_')),'.pdf']
    %Fig = figure('PaperSize',[10*length(Nicknames) Height],'resize','off', 'visible','on','Units','points');
        Fig = figure('PaperUnits','inches','PaperSize',[3*length(Nicknames) 6],'Units','points','resize','on', 'visible','on','DefaultAxesFontSize', 4);
set(0,'defaultAxesFontSize',4)
      set(0, 'DefaultFigureRenderer', 'painters');
          set(gcf,'defaultLegendAutoUpdate','off')
          set(Fig,'defaultAxesColorOrder',Palette)
%   set(Fig,'defaultAxesFontSize',6)
 %     set(Fig, 'DefaultFigureRenderer', 'painters');


% NUMBER CELLS ON
% subplot(2,5,1); 
%plotBoxplot(NumOn,Nicknames,ExpLabels,Jitter,BarW,'# active cells',FontSize,5,Palette,FaceAlpha,LineWidth,FontSizeTitle,[0,100])




% PERC CELLS ON
subplot(2,5,1); 
plotBoxplot(PercOn*100,Nicknames,ExpLabels,Jitter,BarW,'% active cells',FontSize,DotSize,Palette,FaceAlpha,LineWidth,FontSizeTitle,[0,100])  

% PROPORTION SUSTAINED
subplot(2,5,2); 
plotBoxplot(PercCont*100,Nicknames,ExpLabels,Jitter,BarW,'% sustained profile',FontSize,DotSize,Palette,FaceAlpha,LineWidth,FontSizeTitle,[0,100])

% AMPLITUDE
subplot(2,5,3); 
plotViolin(BAmpAll,Nicknames,ExpLabels,100,BarW,'burst amplitude (AU)',FontSize,FaceAlpha,LineWidth,FontSizeTitle,100,Log)

% TIME OFF
subplot(2,5,4); 
plotViolin(OffTimeAll,Nicknames,ExpLabels,0.5,BarW,'off time (min)',FontSize,FaceAlpha,LineWidth,FontSizeTitle,1,Log)

% LENGTH
subplot(2,5,5); 
plotViolin(BLenAll,Nicknames,ExpLabels,0.5,BarW,'burst length (min)',FontSize,FaceAlpha,LineWidth,FontSizeTitle,1,Log)

% ONSET/ENDS VIOLIN
subplot(2,5,6); 
plotViolin(OnsetAll,Nicknames,ExpLabels,2,BarW,'onset/ends (min)',FontSize,0.2,LineWidth,FontSizeTitle,1,0)
plotViolin(EndAll,Nicknames,ExpLabels,2,BarW,'onset/ends (min)',FontSize,FaceAlpha,LineWidth,FontSizeTitle,1,0)

% NUMBER OF BURSTS
subplot(2,5,7); 
plotViolin(BNumAll,Nicknames,ExpLabels,1,BarW,'bursts # / cell',FontSize,FaceAlpha,LineWidth,FontSizeTitle,1,0)

% ONSET/ENDS VIOLIN
subplot(2,5,8); 
plotViolin(TmAll,Nicknames,ExpLabels,100,BarW,'total mRNA (AU)',FontSize,FaceAlpha,LineWidth,FontSizeTitle,1,0)

% MAX LENGTH EACH CELL
subplot(2,5,9); 
plotViolin(MaxLenAll,Nicknames,ExpLabels,0.5,BarW,'max burst length (min)',FontSize,FaceAlpha,LineWidth,FontSizeTitle,1,Log)

% PERIOD
% subplot(2,5,9); 
% plotViolin(BPerAll,Nicknames,ExpLabels,0.5,BarW,'burst period (min)',FontSize,FaceAlpha,LineWidth,FontSizeTitle,1,Log)

% SIZE
subplot(2,5,10); 
plotViolin(BSizeAll,Nicknames,ExpLabels,100,BarW,'burst size (AU)',FontSize,FaceAlpha,LineWidth,FontSizeTitle,1,Log)


%export_fig FileOut -pdf
%saveas(Fig,FileOut,'pdf')
print(Fig,FileOut,'-fillpage','-dpdf');
%print(Fig,FileOut,'-dtiff');
close all
%end

%%

%asterisk(p) returns * ** or *** for each p value depending on
%predetermined levels of significance, can be changed inside asterisk
%function

% a test decision for the null hypothesis that the data in vectors x and y
% comes from independent random samples from normal distributions with
% equal means and equal but unknown variances, using the two-sample t-test.
% The alternative hypothesis is that the data in x and y comes from
% populations with unequal means. The result h is 1 if the test rejects the
% null hypothesis at the 5% significance level, and 0 otherwise
% 1 if coming from different normal distributions
%[h,p,ci,stats] = ttest2(BAmpAll(:,1),BAmpAll(:,2),'Vartype','unequal'); disp(h);p
%[h,p,ci,stats] = ttest2(BAmpAll(:,2),BAmpAll(:,3),'Vartype','unequal'); disp(h);p


% a test decision for the null hypothesis that the data in vector x comes
% from a distribution in the normal family, against the alternative that it
% does not come from such a distribution, using a Lilliefors test.
% 1 if not normally distributed
% lillietest(BAmpAll(:,1))
% lillietest(BAmpAll(:,2))
% lillietest(BAmpAll(:,3))


%a test decision for the null hypothesis that the data in vectors x1 and x2
%are from the same continuous distribution, using the two-sample
%Kolmogorov-Smirnov test. The alternative hypothesis is that x1 and x2 are
%from different continuous distributions. The result h is 1 if the test
%rejects the null hypothesis at the 5% significance level, and 0 otherwise
fileID = fopen([PathMain,'/comp_bursting_violin',cell2mat(join(unique(Selections),'_')),'_',cell2mat(join(Nicknames,'_vs_')),'_ks.txt'],'a');
fprintf(fileID, '\nupdated on \t%s\n',date);
[h,p,ks2stat] = kstest2(BAmpAll(:,1),BAmpAll(:,2));
fprintf(fileID, 'amp12\t%E\t%s\n',p,asterisk(p));
[h,p,ks2stat] = kstest2(OffTimeAll(:,1),OffTimeAll(:,2));
fprintf(fileID, 'off12\t%E\t%s\n',p,asterisk(p));
[h,p,ks2stat] = kstest2(BLenAll(:,1),BLenAll(:,2));
fprintf(fileID, 'length12\t%E\t%s\n',p,asterisk(p));
[h,p,ks2stat] = kstest2(BSizeAll(:,1),BSizeAll(:,2));
fprintf(fileID, 'size12\t%E\t%s\n',p,asterisk(p));
if length(Nicknames) == 3
    [h,p,ks2stat] = kstest2(BAmpAll(:,2),BAmpAll(:,3));
    fprintf(fileID, 'amp23\t%E\t%s\n',p,asterisk(p));
    [h,p,ks2stat] = kstest2(OffTimeAll(:,2),OffTimeAll(:,3));
    fprintf(fileID, 'off23\t%E\t%s\n',p,asterisk(p));
    [h,p,ks2stat] = kstest2(BLenAll(:,2),BLenAll(:,3));
    fprintf(fileID, 'length23\t%E\t%s\n',p,asterisk(p));
    [h,p,ks2stat] = kstest2(BSizeAll(:,2),BSizeAll(:,3));
    fprintf(fileID, 'size23\t%E\t%s\n',p,asterisk(p));
    [h,p,ks2stat] = kstest2(BAmpAll(:,1),BAmpAll(:,3));
    fprintf(fileID, 'amp13\t%E\t%s\n',p,asterisk(p));
    [h,p,ks2stat] = kstest2(OffTimeAll(:,1),OffTimeAll(:,3));
    fprintf(fileID, 'off13\t%E\t%s\n',p,asterisk(p));
    [h,p,ks2stat] = kstest2(BLenAll(:,1),BLenAll(:,3));
    fprintf(fileID, 'length13\t%E\t%s\n',p,asterisk(p));
    [h,p,ks2stat] = kstest2(BSizeAll(:,1),BSizeAll(:,3));
    fprintf(fileID, 'size13\t%E\t%s\n',p,asterisk(p));
end
if length(Nicknames) == 4
    [h,p,ks2stat] = kstest2(BAmpAll(:,3),BAmpAll(:,4));
    fprintf(fileID, 'amp34\t%E\t%s\n',p,asterisk(p));
    [h,p,ks2stat] = kstest2(OffTimeAll(:,3),OffTimeAll(:,4));
    fprintf(fileID, 'off34\t%E\t%s\n',p,asterisk(p));
    [h,p,ks2stat] = kstest2(BLenAll(:,3),BLenAll(:,4));
    fprintf(fileID, 'length34\t%E\t%s\n',p,asterisk(p));
    [h,p,ks2stat] = kstest2(BSizeAll(:,3),BSizeAll(:,4));
    fprintf(fileID, 'size34\t%E\t%s\n',p,asterisk(p));
end

fclose(fileID);
end