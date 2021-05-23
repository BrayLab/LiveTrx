function[Baseline, baseline,MaxFBG, BG, BGLabels] = BaselinefromBG(MaxF,FoldOverBaseline,Labels, Smooth,minOn)
    %% define baseline
    %Fig1 = figure('PaperSize',[50 12],'PaperUnits','inches','resize','on', 'visible','on');
    %Fig1.Renderer='Painters';
    MeanBaseline = nanmean(MaxF,2);
    %plot(MeanBaseline); hold on
    LFit = fitlm([1:length(MeanBaseline)],MeanBaseline);
    baseline = LFit.Coefficients.Estimate'
    Baseline = [1:length(MeanBaseline)]*baseline(2) + baseline(1);
    %figure(Fig1);subplot(131);plot(TimeScale,MeanBaseline);hold on; plot (TimeScale,Baseline);ylim([0,2000])
     %   xlabel('time into nc14 (min)'); ylabel('Mean Fluorescence (AU)'); legend('show'); legend boxoff; box off
        % fit to line
    % save baseline
    %% get background and redefine baseline
    %Smooth = 3;
    MedFilt = medfilt1(MaxF,Smooth,[],1,'includenan','zeropad');
    MedFilt(MedFilt==0) = NaN;
    %plot(MedFilt); hold on
    %plot (Baseline);
    % onoff
    %minOn = 5;
    OnOff = MedFilt > [Baseline'.*FoldOverBaseline];
    %plot(MaxF.*OnOff);

    % get background
    %OnOffString = join(string(double([(~isnan(MaxF))])),'',1);
    OnOffString = join(string(double(OnOff)),'',1);
    BG = cellfun(@(x) isempty(x), strfind(OnOffString,repmat('1',1,minOn)));
    BGLabels = Labels (BG);
   
    %%
    MaxFBG = MaxF(:,BG);

    if FoldOverBaseline ~=0
        % redefine baseline
         %figure(Fig1);subplot(132);plot(TimeScale, nanmean(MaxFBG,2)); hold on; plot(TimeScale,nanmean(MaxF(:,~BG),2));ylim([0,2000])
          %   xlabel('time into nc14 (min)'); ylabel('Mean Fluorescence (AU)'); legend('show'); legend boxoff; box off
        MeanBaseline = nanmean(MaxFBG,2);
        %plot(MeanBaseline); hold on
        LFit = fitlm([1:length(MeanBaseline)],MeanBaseline);
        baseline = LFit.Coefficients.Estimate'
        Baseline = [1:length(MeanBaseline)]*baseline(2) + baseline(1);
        %figure(Fig1);subplot(133);plot(TimeScale,MeanBaseline);hold on; plot (TimeScale,Baseline);ylim([0,2000]); 
        %xlabel('time into nc14 (min)'); ylabel('Mean Fluorescence (AU)'); legend('show'); legend boxoff; box off
        %print(Fig1,['~/Google Drive jf565/BG.pdf'],'-fillpage', '-dpdf');

        %plot (Baseline);
    end
end