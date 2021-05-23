function[MedFilt, OnOff,Properties, Fig] = DefineExpAll(MaxF,MaxFBG,CentX, Labels, Baseline, FoldOverBaseline,TimeRes,Smooth,SplitShortLong,SplitEarly,nc14,Delay,minOn)
%% Define expression all
% onoff
MedFilt = medfilt1(MaxF,Smooth,[],1,'includenan','zeropad');
MedFilt(MedFilt==0) = NaN;
OnOff = MedFilt > (Baseline'*FoldOverBaseline);
OnOff = CleanOnOff(OnOff,minOn);
[OnOff] = CleanNaNs(MedFilt,OnOff, minOn*2);
%plot(MaxF.*OnOff);
% AP position, mean all CentX
AP = nanmean(CentX,1);
Properties = table(Labels, AP','VariableNames',{'Label','AP_position'});
% first and last track
% early and late exp
% properties: on off, length on, type
SplitMSE = max(round(nc14-Delay + SplitEarly*60/TimeRes),1); %frame to split
%SplitLate = round(nc14-Delay + SplitShortLong*60/TimeRes); %frame to split
OnOffString = join(string(double(OnOff(SplitMSE:end,:))),'',1);
%OnOffString = join(string(double([(~isnan(MaxF(SplitMSE:end,:)))])),'',1);
OnPeriods = strfind(OnOffString,repmat('1',1,minOn));
Midline = cellfun(@(x) ~isempty(x), OnPeriods);
MidlineLabels = Labels(Midline);
On = cellfun(@(x) x(1), OnPeriods(Midline));
Off = cellfun(@(x) x(end), OnPeriods(Midline))+minOn-1;
Properties.Onset = nan(size(Properties.Label)); 
Properties.End = nan(size(Properties.Label));
Properties.TimeOn = nan(size(Properties.Label));
Properties.Onset(Midline) = [(On + SplitMSE -nc14 + Delay).*TimeRes./60]';
Properties.End(Midline) =[(Off + SplitMSE -nc14 + Delay).*TimeRes./60]';
Properties.TimeOn = Properties.End - Properties.Onset;
OnOffStringE = join(string(double([(~isnan(MaxF(1:SplitMSE,:)))])),'',1);
OnEarly = ~cellfun(@(x) isempty(x), strfind(OnOffStringE,repmat('1',1,minOn)));
EarlyLabels = Labels (OnEarly);
EarlyOnly = OnEarly & ~Midline; EarlyOnlyLabels = Labels (EarlyOnly);
%% split short/long and write type of expression
LongMid = [Properties.End > SplitShortLong]; LongMidLabels = Labels(LongMid);
ShortMid = [Properties.End < SplitShortLong]; ShortMidLabels = Labels(ShortMid);
Properties.Type = string(repmat('NaN',length(Properties.Label),1));
%Properties.Type(BG) = string(repmat('BG',length(BGLabels),1));
Properties.Type(EarlyOnly) = string(repmat('EarlyOnly',length(EarlyOnlyLabels),1));
Properties.Type(ShortMid) = string(repmat('ShortMidline',length(ShortMidLabels),1));
Properties.Type(LongMid) = string(repmat('LongMidline',length(LongMidLabels),1));
% add late track and gastrulated
% late track
NaNString = join(string(double(~isnan(MaxF))),'',1);
    Beginings = strfind(NaNString,'1');
    StartF = cellfun(@(x) x(1)+1, Beginings);
    LateTrack = (StartF-nc14+Delay).*TimeRes./60 > 40; LateTrackLabels = Labels(LateTrack);
   Properties.Type(LateTrack) = string(repmat('LateTrack',length(LateTrackLabels),1));
% gastrulated

%% split exp
MMidline = nanmean(MaxF(:,Midline),2);
MShortMid = nanmean(MaxF(:,ShortMid),2);
MLongMid = nanmean(MaxF(:,LongMid),2);
MBG = nanmean(MaxFBG,2);
MEarly = nanmean(MaxF(:,EarlyOnly),2);
TimeScale = ([1:length(Baseline)]-nc14+Delay)*TimeRes/60;
Fig = figure('Visible','off');
plot(TimeScale,MMidline); hold on;
plot(TimeScale,MShortMid); hold on;
plot(TimeScale,MLongMid); hold on;
plot(TimeScale,MBG);
plot(TimeScale,MEarly);
plot(TimeScale,Baseline);
xlabel('Time into nc14 (min)');
ylabel('Mean F (AU)');
end