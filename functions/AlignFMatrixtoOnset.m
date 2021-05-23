function[NormAligned,TimeScaleAligned] = AlignFMatrixtoOnset(NormMerged,PropertiesMerged,TimeScaleMerged,MaxTime,TimeRes)
    OnsetsF = [PropertiesMerged.Onset]*60/TimeRes - 3;
    NormAligned = nan(round(MaxTime*60/TimeRes),size(NormMerged,2));
    DelayMerged = find(round(TimeScaleMerged,4) == 0);
    %if isempty(DelayMerged); DelayMerged = 0; end
    for l = 1:length(OnsetsF)
        if isnan(OnsetsF(l)) == 0
            Start = round(OnsetsF(l) + DelayMerged);
            End = min(MaxTime*60/TimeRes+Start-1,size(NormMerged,1));
            NormAligned(1:End-Start+1,l) = NormMerged(Start:End,l);
        end
    end
    %figure; plot(nanmean(NormAligned,2))
    TimeScaleAligned = [1:size(NormAligned,1)].*TimeRes./60;
    
end