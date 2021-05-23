function[] = CompareMembranes(SelectedN,SelectedNames,Which,Info,Path,Exps,Channels, Palette, Boo)

% COMPARE MEMBRANES


Growth = Boo(1); BoxplotsLength = Boo(2); CorrCell = Boo(3);


mkdir(Path)

if BoxplotsLength | CorrCell
    
    %Palette2 = kron(Palette,ones(2,1)); %duplicate each row consecutively
    set(0,'defaultAxesColorOrder',Palette)
    set(0,'defaultLegendAutoUpdate','on')
    
    ToSave = [Path,char(join(SelectedN,'_vs_'))];

    if CorrCell
        Fig = figure('PaperSize',[25 20],'PaperUnits','inches','resize','on', 'visible','on');
    end
    
    MeanOn = cell(1,length(SelectedN));
    MembTime1 = cell(1,length(SelectedN));
    MembTime2 = cell(1,length(SelectedN));
    MembTime3 = cell(1,length(SelectedN));
    MembTime4 = cell(1,length(SelectedN));
    MembTime5 = cell(1,length(SelectedN));


    for i = 1:length(SelectedN)
        Selection = '';
            %try
                Index = find(cellfun(@(x) strcmp(x,SelectedN{i}),table2array(Exps(:,Which)))==1)';
                for x = 1:length(Index)
                    %try
                    Experiment = [Exps.Nickname{Index(x)},' ',num2str(Exps.Rep(Index(x)))]
                    PathToSave = [Info.Path{Index(x)},Info.File{Index(x)},...
                    Info.Name{Index(x)},Info.File{Index(x)}];                 
                    load([PathToSave,'_Data.mat']);
                    Properties = Data.Properties;
                    MeanOn{i} = [MeanOn{i}, quantile(Properties.Onset,0.25)];
                    %MeanOn{i} = [MeanOn{i}, NaN];
                    %MeanOn(x) = mean(Properties.Onset, 'omitnan');
                    MembData = readtable([PathToSave,'_membranes.txt']);
                    MembTime1{i} = [MembTime1{i},MembData.Time(1)];   
                    MembTime2{i} = [MembTime2{i},MembData.Time(2)];   
                    MembTime3{i} = [MembTime3{i},MembData.Time(3)];   
                    MembTime4{i} = [MembTime4{i},MembData.Time(4)];   
                    MembTime5{i} = [MembTime5{i},MembData.Time(5)];             
                    %end
                end 
           % end
    if CorrCell
        subplot(221); hold on
        PlotScatter(MembTime2{i},MeanOn{i},SelectedNames{i},'Cellularization Time (min)', 'Q1 Onset (min)',1,0); 

        subplot(222); hold on
        PlotScatter(MembTime3{i},MeanOn{i},SelectedNames{i},'Cellularization Time (min)', 'Q1 Onset (min)',1,0); 
        
        subplot(223); hold on
        PlotScatter(MembTime4{i},MeanOn{i},SelectedNames{i},'Cellularization Time (min)', 'Q1 Onset (min)',1,0); 

        subplot(224); hold on
        PlotScatter(MembTime5{i},MeanOn{i},SelectedNames{i},'Cellularization Time (min)', 'Q1 Onset (min)',1,0); 
    end

    end

    MeanOnAll = Cell2Mat(MeanOn);
    MembTime1All = Cell2Mat(MembTime1);
    MembTime2All = Cell2Mat(MembTime2);
    MembTime3All = Cell2Mat(MembTime3);
    MembTime4All = Cell2Mat(MembTime4);
    MembTime5All = Cell2Mat(MembTime5);

    if CorrCell
                subplot(221); hold on
        PlotScatter(MembTime2All,MeanOnAll,'All','Cellularization Time (min)', 'Q1 Onset (min)',0,1); 
        xlim([5,55])
        ylim([5,55])
                subplot(222); hold on
        PlotScatter(MembTime3All,MeanOnAll,'All','Cellularization Time (min)', 'Q1 Onset (min)',0,1); 
        xlim([5,55])
        ylim([5,55])
        
        subplot(223); hold on
        PlotScatter(MembTime4All,MeanOnAll,'All','Cellularization Time (min)', 'Q1 Onset (min)',0,1); 
        xlim([5,55])
        ylim([5,55])

        subplot(224); hold on
        PlotScatter(MembTime5All,MeanOnAll,'All','Cellularization Time (min)', 'Q1 Onset (min)',0,1); 
        xlim([5,55])
        ylim([  5,55])

        print(Fig,[ToSave,'_MembCorr.pdf'],'-fillpage', '-dpdf');
        close all
    end

    % BOXPLOT CELLULARIZING TIMES
    if BoxplotsLength
        Fig = figure('PaperSize',[30 15],'PaperUnits','inches','resize','on', 'visible','on');

        SelectedNamesX = repmat(SelectedNames(:),1,5);
        SelectedNamesX = SelectedNamesX(:);
        Palette4 = [Palette;Palette;Palette;Palette;Palette];
        set(0,'defaultAxesColorOrder',Palette)

        Jitter = 0.5; %Jitter./2 cant be > BarW 
        BarW = 0.35; FaceAlpha = 0.3; DotSize = 15; LineWidth = 2; FontSizeTitle = 16; FontSize = 14;
        set(gca,'FontSize',FontSize)
         %subplot(1,10,[6:7]); hold on     
        plotBoxplot([MembTime1All,MembTime2All,MembTime3All,MembTime4All,MembTime5All],SelectedNamesX,SelectedNamesX,Jitter,BarW,'',FontSize,DotSize,Palette4,FaceAlpha,LineWidth,FontSizeTitle,[0,60]);
         ylabel('Cellularization time (min)')

        print(Fig,[ToSave,'_MembTimesBoxplot.pdf'],'-fillpage', '-dpdf');
        close all
    end
end


%%

% Which = 1 % 1 to merge all repeats, 3 for each
% UniqueN = table2array(unique(Exps(:,Which))); 
% Subplots = 14;
% XLim = [0,60]
% for n = 1:(floor(length(UniqueN)/Subplots)+1)
%     Fig1 = figure('PaperSize',[50 50],'PaperUnits','inches','resize','on', 'visible','on');
%     Fig1.Renderer='Painters';
%     %
%     Selection = 'MSE';
%     %try
%     ToSave = [Path];
%     for i = 1:Subplots
%         try
%             Index = find(cellfun(@(x) strcmp(x,UniqueN{(n-1)*Subplots+i}),table2array(Exps(:,Which)))==1)';
%             for x = 1:length(Index)
%                 
%                 Experiment = [Exps.Nickname{Index(x)},' ',num2str(Exps.Rep(Index(x)))]
%                 PathToSave = [Info.Path{Index(x)},Info.File{Index(x)},...
%                     Info.Name{Index(x)},Info.File{Index(x)}];                 
%                 load([PathToSave,'_Data.mat']);
%                 MaxF = Data.MaxF;
%                 CentY = Data.CentY;                
%                 CentX = Data.CentX;
%                 CentZ = Data.CentZ;
%                 %OnOff = Data.OnOff;
%                 Properties = Data.Properties;
%                 Baseline = Data.Baseline;
%                 TimeScale = Data.TimeScale;
%                 Bits = Info.Bits(Index(x));
%                 XYRes = Info.XYRes(Index(x));
%                 ZRes = Info.ZRes(Index(x));
%                 TimeRes = Info.TimeRes(Index(x));
%                 Norm = (MaxF-Baseline').*Baseline(1)./Baseline';
%                 Norm(Norm==0) = NaN;
%                 if strcmp(Selection,'') == 1
%                     Selected = [Properties.Type ~= 'BG'];
%                 else
%                     Selected = [Properties.Type ~= 'EarlyOnly' & Properties.Region == Selection];
%                 end
%                 TimeRes = Info.TimeRes(Index(x));
%                 nc14 = Info.nc14(Index(x));
%                 Delay = Info.Delay(Index(x));
%                 MovX = abs(diff(CentX,1,1))*XYRes;
%                 MovY = abs(diff(CentY,1,1))*XYRes;
%                 MovZ = abs(diff(CentZ,1,1))*ZRes;
%                 MovXYZ = sqrt(MovX.^2+MovY.^2+MovZ.^2)./TimeRes;
%                 figure(Fig1);
%                 set(Fig1,'defaultAxesColorOrder',[[144,191,91]./255;[209,28,71]./255]);
%                 subplot(Subplots/2,2,i); hold on
%                 yyaxis('left')
%                 plot(TimeScale,nanmean(Norm,2),'-', 'Color',[144,191,91]./255,'LineWidth',1,'LineStyle','-')
%                 %plot(TimeScale,nanmean(-CentY,2),'LineStyle','-')
%                 %plot(TimeScale,nanmean(-CentY,2)+nanstd(CentY,1,2),'LineStyle','--')
%                 %plot(TimeScale,nanmean(-CentY,2)-nanstd(CentY,1,2),'LineStyle','--')
%                 ylim([-300,1500])
%                 ylabel('Mean F (AU)')
%                 yyaxis('right')
%                 ylim([0,30]) 
%                 try
%                     MembData = readtable([PathToSave,'_membranes.txt']);
%                     MembLength = MembData.MembLenght*50*XYRes;
%                     plot(MembData.Time, MembLength,'.-','Color',[209,28,71]./255,'LineWidth',1, 'LineStyle','-')
%                     %stairs(MembData.Time(1:end-1), diff(MembLength)./(diff(MembData.Time)*60))
%                 end
%                 xlim([XLim]); 
%                 xlabel('Time into nc14 (min)')
%                 ylabel('Membrane length (um)')
%                 title(Experiment)
%                 
%                
%             end 
%     end
% 
%     end
%     if n==1
%           print(Fig1,[ToSave,Selection,MetaFile,'memb.ps'],'-fillpage', '-dpsc');
%      else
%           print(Fig1,[ToSave,Selection,MetaFile,'memb.ps'],'-fillpage', '-dpsc','-append');
%      end
% end
% close all

%%
%%  COMPARE LENGTH OF MEMBRANES

% 
% MetaFile = '';
% MetaFile = ' membranes';
% 
% 
% Path = ['/Users/julia/Google Drive jf565/comp MS2/',MetaFile,'/'];
% mkdir(Path)
% 
% 
% Info = readtable(['~/Google Drive jf565/MATLAB_R_scripts/metadata MS2 3D',MetaFile,'.txt'],'ReadVariableNames', true,'Delimiter','\t');
% Exps = cell2table(Info.Nickname, 'VariableNames',{'Nickname'});
% Exps.Rep = Info.Rep;
% for x=1:size(Exps,1); Exp{x} = [Exps.Nickname{x},' ',num2str(Exps.Rep(x))];end
% Exps.Exp = Exp';
% % RGB values from ggplot colors
% %PalettePoints <-c('steelblue1','olivedrab2',"orangered2",'darkgoldenrod1','grey70','grey35','slateblue3')
% %PaletteMeans <-c('steelblue4','olivedrab4',"orangered3",'darkgoldenrod3','grey40','grey20','slateblue4')
% PalettePoints = [179,238,58;99,184,255;238,64,0;255,185,15;179,179,179;89,89,89;105,89,205;100,100,100]./255;
% PaletteMeans = [105,139,34;54,100,139;205,55,0;205,149,12;102,102,102;51,51,51;71,60,139;80,80,80]./255;
% 
% %PalettePoints = [142,183,36;61,131,183;105,139,34;54,100,139;62,81,16;33,63,86]./255;
% %PaletteMeans = [142,183,36;61,131,183;105,139,34;54,100,139;62,81,16;33,63,86]./255;
% 
% 
% set(0,'defaultAxesFontSize',16)
% %set(0,'defaultAxesFontName', 'Myriad Pro')
% 
% %
% 
% Which = 1 % 1 to merge all repeats, 3 for each
% UniqueN = table2array(unique(Exps(:,Which))); 
% 
% PaletteMain = [130,130,130;80,80,80; %grey1, grey2. 1,2
%     142,183,36;105,139,34;62,81,16; %green1,green2,green3. 3:5
%     61,131,183;54,100,139;33,63,86; %blue1,blue2,blue3. 6:8
%     250,174,64;226,140,26; %yellow1, yellow2. 9:10
%     167,51,170;104,15,107; %pink1, pink2. 11:12
%     78,69,139; 54, 44,96; %purple1, purple2. 13:14
%     205,55,0]./255; %red. 15
% 
% PairstoSelect = {
%     {'m5m8peve gap43mCh','spiderGFP hisRFP','NiGFP','m5m8peve DlmSca'},...
%         {'shgGFP DlmSca','shgGFP DlmSca'},...
%         {'m5m8peve gap43mCh','spiderGFP hisRFP'},...
%         {'m5m8peve DlmSca','slam:+ DlmSca','slam DlmSca'}
%         }
%     
% Nicknames = {
%     {'gap43','spider','N','Dl'},...
%     {'shg','Dl',},...
%     {'gap43','spider'},...
%     {'+','slam^{+/+/-}','slam^{-/-}'},...
%         }
% PaletteEach = {[1,2,3,11],[3,11],[9,10],[1,12,11]}
% 
% ChannelsEach = {[2,1,1,2],[1,2],[2,1],[2,1,1]}
% 
% n = 1

if Growth

    %SelectedN = PairstoSelect{n}
    %SelectedNames = Nicknames{n}
    %Channels = ChannelsEach{n}
    Palette2 = kron(Palette,ones(2,1)); %duplicate each row consecutively
    set(0,'defaultAxesColorOrder',Palette2)
    set(0,'defaultLegendAutoUpdate','on')


    ToSave = [Path,char(join(SelectedN,'_vs_'))];
            Fig1 = figure('PaperSize',[20 50],'PaperUnits','inches','resize','on', 'visible','on');
            %Fig2 = figure('PaperSize',[20 50],'PaperUnits','inches','resize','on', 'visible','on');

    XLim = [-20,50];
    YLimits = [0,40];

    for i = 1:length(SelectedN)
        Channel = str2num(Channels{i});
        MembDataMerged = [];
        TimeScaleMerged = [];
        PropertiesMerged = table();
        PropertiesMerged.Type = zeros(0);
        Selection = '';
           % try
                Index = find(cellfun(@(x) strcmp(x,SelectedN{i}),table2array(Exps(:,Which)))==1)';
                for x = 1:length(Index)
                    %try
                    Experiment = [Exps.Nickname{Index(x)},' ',num2str(Exps.Rep(Index(x)))]
                    PathToSave = [Info.Path{Index(x)},Info.File{Index(x)}]; 
                    Parameters = Info(Index(x),:);
                    Table2Vars(Parameters);
                    TimeScale = ([1: Frames] - nc14+Delay).*TimeRes./60;
                    %MeanOn(x) = mean(Properties.Onset, 'omitnan');
                    if Channel == 1 & strcmp(Channel1,'memb')                    
                        MembData = readtable([PathToSave,'_MembraneLengthCh1.txt']);         
                    end
                    if Channel == 2 & strcmp(Channel2,'memb')                    
                        MembData = readtable([PathToSave,'_MembraneLengthCh2.txt']);          
                    end
                    MembData = MembData.Var1.*XYRes;
                    Properties = table(NaN);
                    Properties.Type = 'ShortMidline';
                    [MembDataMerged,PropertiesMerged,TimeScaleMerged] = MergeFMatrix(MembData,MembDataMerged,Properties,PropertiesMerged,TimeScale,TimeScaleMerged,TimeRes);
                    
                end 
            %end
            
           ColorArg = [Palette(i,:);Palette(i,:);Palette(i,:);Palette(i,:);Palette(i,:);Palette(i,:);Palette(i,:);Palette(i,:)];
                %plot(TimeScaleMerged,MembDataMerged); hold on
                Fig1 = PlotMeansFractionShaded(MembDataMerged,MembDataMerged,MembDataMerged,TimeScaleMerged,ones(size(MembDataMerged,2),1),PropertiesMerged,Bits,ColorArg,Fig1,1,SelectedN{i},XLim, YLimits);
                figure(Fig1)
                subplot(411); ylabel('Membrane length (um)')
                subplot(412); ylabel('Membrane length (um)')
                subplot(413); ylabel('Membrane length (um)')
                subplot(414); ylabel('Membrane length (um)')

    %             SpeedMemb = diff(MembDataMerged,1,1)./TimeRes .*60;
    %             Fig2 = PlotMeansFractionShaded(SpeedMemb,SpeedMemb,SpeedMemb,TimeScaleMerged,ones(size(SpeedMemb,2),1),PropertiesMerged,Bits,ColorArg,Fig2,1,SelectedN{i},XLim, YLimits);
    %             figure(Fig2)
    %             subplot(411); ylabel('Membrane growth rate (um/min)')
    %             subplot(412); ylabel('Membrane growth rate (um/min)')
    %             subplot(413); ylabel('Membrane growth rate (um/min)')
                subplot(414); ylabel('Membrane growth rate (um/min)')


    end

            print(Fig1,[ToSave,Selection,'_means.pdf'],'-fillpage', '-dpdf');
            %print(Fig2,[ToSave,Selection,'_means.pdf'],'-fillpage', '-dpdf');
end
%close all