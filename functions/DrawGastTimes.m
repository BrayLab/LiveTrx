function[Data] = DrawGastTimes(Data,PathData,Parameters,Selection,Confirm,varargin)
%will always ask to draw if it doesnt exist, or not in right region
% but will only show plot to confirm if Confirm = 1

        Table2Vars(Parameters);
        %PathToSave = [Path,File,Name,File]; 
        load(PathData);
        Struct2Vars(Data);
                
        try
        if ~isempty(varargin{1})
            SplitEarly = varargin{1}; 
        end
        end
            SelectionToSave = Selection;
            if strcmp(Selection,'') == 1
                Selected = [Properties.Type ~= 'BG'];
            else
                Selected = Properties.Type ~= 'EarlyOnly' & Properties.Region == Selection;
                if contains(Selection,'|')
                    Selected = [];
                    Selections = strsplit(Selection,'|');
                    SelectionToSave = join(Selections,'');SelectionToSave = SelectionToSave{:};
                    for s = 1:length(Selections)
                    	Selected(s,:) = Properties.Type ~= 'EarlyOnly' & Properties.Region == Selections{s};
                    end
                    Selected = sum(Selected,1) ~= 0;      
                end    
            end
                
                minOn = 5;
                Labels = Properties.NewLabel;
                FoldOverBaseline = 1.2;
                [~, ~,PropertiesNew, ~] = DefineExpAll(MaxF,MaxFBG,CentX,Labels,Baseline,FoldOverBaseline,TimeRes,3,60,SplitEarly,nc14,Delay,minOn);
                %close Fig
                Properties.Onset = PropertiesNew.Onset;
                Properties.Type = PropertiesNew.Type;
                Norm = ((MaxF-Baseline').*Baseline(1)./Baseline').*(2.^(12-Bits));
                Norm(Norm==0) = NaN;
                NormBG = (MaxFBG-Baseline').*Baseline(1)./Baseline';
                NormBG(NormBG==0) = NaN;
%                MeantoPlot = [nanmean(Norm(1:max(1,nc14-Delay+15*60/TimeRes),Properties.Type == 'EarlyOnly'),2);...
                   % nanmean(Norm(max(1,nc14-Delay+15*60/TimeRes)+1:end,Selected),2)];
                MeantoPlot = nanmean(Norm(:,Selected),2);

                %YPosition = nanmean(CentY(:,Selected),2);
                %MSEPos = nanmean(YPosition(max(1,nc14-Delay+15*60/TimeRes):max(1,nc14-Delay+30*60/TimeRes)));
                %YPosition = YPosition-MSEPos;
                MovX = abs(diff(CentX,1,1))*XYRes;
                MovY = abs(diff(CentY,1,1))*XYRes;
                MovZ = abs(diff(CentZ,1,1))*ZRes;
                MovXYZ = sqrt(MovX.^2+MovY.^2+MovZ.^2);
                Mov = MovXYZ./TimeRes;
                Data.MovY = MovY;
                Data.MovX = MovX;
                Data.MovZ = MovZ;
                Data.MovXYZ = MovXYZ;

                MovSelected = Mov(:,Selected);
                
                GastValues = [];
                Answer = 'No';
                pause(0.5)
                try
                    GastValues = Data.(['GastValues',Selection]);
                    PeakT = GastValues(:,1);
                    PeakF = GastValues(:,2);
                    GastT = GastValues(:,3);
                    
                    if Confirm == 1
                        FigDraw = figure();
                        subplot(211)
                        plot(TimeScale,MeantoPlot, '-','Color',[144,191,91]./255,'LineWidth',1); hold on 
                        title('Select Peak#1, transition and Peak#2 (NaN @ t<0)')
                        %[PeakT,PeakF] = ginput(3); 
                        %PeakF(PeakT<0) = NaN; PeakT(PeakT<0) = NaN;
                        plot(PeakT,PeakF,'*r')

                        subplot(212)
                        yyaxis left
                        plot(TimeScale,CentY(:,Selected), '-','Color',[144,191,91]./255,'LineWidth',1); hold on 
                        %plot(TimeScale,YPosition, '-','Color',[144,191,91]./255,'LineWidth',2); hold on 
                        title('Select Begining, peak and end of gastrulation (NaN @ t<0)')
                        yyaxis right
                        plot(TimeScale(1:end-1),nanmean(MovSelected,2),'-','Color',[209,28,71]./255,'LineWidth',1)
                        %plot(TimeScale(1:end-1),medfilt1(nanmean(MovSelected,2),5),'-','Color',[209,28,71]./255,'LineWidth',1)
                        %[GastT,~] = ginput(3); 
                        %GastT(GastT<0) = NaN;
                        yyaxis left
                        plot([GastT,GastT],[0,400],'--')
                        hold off
                        Answer = questdlg('OK?');
                    else
                        Answer = 'Yes';
                    end
                end    
                    
                    % DELETE LATER. JUST TO SAVE PLOT
                    %%%%%%%%%%%%%%%%%%%%%%%%%
%                     Fig2 = figure('PaperSize',[15 15],'PaperUnits','inches','resize','on', 'visible','on');
%                         subplot(211)
%                         plot(TimeScale,MeantoPlot, '-','Color',[144,191,91]./255,'LineWidth',1); hold on 
%                         title('Select Peak#1, transition and Peak#2 (NaN @ t<0)')
%                         plot(PeakT,PeakF,'*r')
%                         xlim([0,70])
% 
%                         subplot(212)
%                         yyaxis left
%                         plot(TimeScale,CentY(:,Selected), '-','Color',[[61,131,183]./255,0.3],'LineWidth',0.5); hold on 
%                         title('Select Begining, peak and end of gastrulation (NaN @ t<0)')
%                         yyaxis right
%                         plot(TimeScale(1:end-1),medfilt1(nanmean(MovXYZ,2),5),'-','Color',[33,63,86]./255,'LineWidth',1)
%                         yyaxis left
%                         plot([GastT,GastT],[0,400],'--')
%                         xlim([0,70])
% 
%                         hold off
                     %%%%%%%%%%%%%%%%%%%%%%%  
%                      mkdir([ToSave,'PlotsEach/'])
%                     print(Fig2,[ToSave,'PlotsEach/',Experiment,'.pdf'],'-fillpage', '-dpdf');
%                 
                
                %catch    
                    %disp('here')

                    while strcmp(Answer,'Yes') ~= 1
                        Fig2 = figure();
                        subplot(211)
                        plot(TimeScale,MeantoPlot, '-','Color',[144,191,91]./255,'LineWidth',1); hold on 
                        xlim([0,inf])
                        title('Select Peak#1, transition and Peak#2 (NaN @ t<0)')
                        [PeakT,PeakF] = ginput(3); 
                        PeakF(PeakT<0) = NaN; PeakT(PeakT<0) = NaN;
                        plot(PeakT,PeakF,'*r')

                        subplot(212)
                        yyaxis left
                        plot(TimeScale,CentY(:,Selected), '-','Color',[144,191,91]./255,'LineWidth',1); hold on 
                        title('Select Begining, peak and end of gastrulation (NaN @ t<0)')
                        yyaxis right
                        %plot(TimeScale(1:end-1),medfilt1(nanmean(MovSelected,2),5),'-','Color',[209,28,71]./255,'LineWidth',1)
                        plot(TimeScale(1:end-1),nanmean(MovSelected,2),'-','Color',[209,28,71]./255,'LineWidth',1)
                        xlim([0,inf])
                        [GastT,~] = ginput(3); 
                        GastT(GastT<0) = NaN;
                        yyaxis left
                        plot([GastT,GastT],[0,400],'--')
                        hold off
                        Answer = questdlg('OK?');
                        if strcmp(Answer,'Cancel' )  
                            close(FigDraw)
                            GastValues = nan(3,3);
                            close all
                            break
                        end
                        %save mat
                        MeanF = [nanmean(MeantoPlot(max(1,30*60/TimeRes+nc14-Delay+1):min(Frames,+50*60/TimeRes+nc14-Delay+1))),NaN,NaN]';
                        GastValues = [PeakT,PeakF,GastT,MeanF];                        
                    end
                    % add to Data
                    Data.(['GastValues',SelectionToSave]) = GastValues;
                    save(PathData,'Data');
                    %pause(0.1)
                %end
  
                %%%%% uncomment to just draw
%                 figure(Fig1);
%                 set(Fig1, 'Visible','off')
%                 subplot(Subplots/NCol,NCol,i); hold on
%                 yyaxis('left')
%                 plot(TimeScale,MeantoPlot, '-','Color',[144,191,91]./255,'LineWidth',1)
%                 plot(TimeScale,nanmean(NormBG,2), '-','Color',[200,200,200]./255,'LineWidth',1)
%                 plot(PeakT,PeakF,'*r')
%                 ylim([-100,1100])
%                 xlim([XLim]); 
%                 ylabel('Mean F (AU)')
%                 
%                 yyaxis('right')
%                 xlim([XLim]); 
%                 ylim([0,0.25]) 
%                 plot(TimeScale(1:end-1),medfilt1(nanmean(MovXYZ(:,Selected),2),5),'-','Color',[209,28,71]./255,'LineWidth',1)
%                 plot([GastT,GastT],[0,0.25],'--','Color',[209,28,71]./255,'LineWidth',1)
%                 title(Experiment)
%                     xlabel('Time into nc14 (min)')
%                     ylabel('Mean DV displacement (um/s)')


%                    
%                 plot(TimeScale,CentY(:,Selected),'LineStyle','-')
%                 plot(TimeScale,nanmean(CentY(:,Selected),2),'LineStyle','-')
%                 plot(TimeScale,nanmean(CentY(:,Selected),2)+nanstd(CentY,1,2),'LineStyle','--')
%                 plot(TimeScale,nanmean(CentY(:,Selected),2)-nanstd(CentY,1,2),'LineStyle','--')
%                 ylim([0,400])

end