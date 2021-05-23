function[BurstNum,BurstLength,BurstPeriod,BurstPeak, BurstMax,OffTime,BurstSize] = CountBursts(MaxF, OnOff, Selected,minOn,StartCount,TimeRes)
    OnOff = CleanOnOff(OnOff(max(StartCount-1,1):end,Selected),minOn);
    % force to start and end in 0 to count all ons and offs. 
    OnOff(1,:) = 0;
    OnOff(end+1,:) = 0;
    MaxF = MaxF(max(StartCount-1,1):end,Selected);
    OnOffString = join(string(double(OnOff)),'',1);
    %OnOffString = join(string(double([(~isnan(MaxF(SplitMSE:end,:)))])),'',1);
    Beginings = strfind(OnOffString,['0',repmat('1',1,minOn)]);
    Ends = strfind(OnOffString,[repmat('1',1,minOn),'0']);
    Beginings = cellfun(@(x) x+1, Beginings,'UniformOutput', false);
    Ends = cellfun(@(x) x+minOn-1, Ends,'UniformOutput', false);
    %OnOffString
    BurstNum = cellfun(@(x) size(x,2), Beginings,'UniformOutput', false);
    
    BurstLength = {}; BurstPeak = {}; BurstPeriod = {}; OffTime = {};BurstMax = {};BurstSize = {};
    % only count burst that have begining and end
    GoodOnes = cellfun(@(x,y) size(x,2)>0, Beginings,Ends);
    %GoodOnes = cellfun(@(x,y) (size(x,2)==size(y,2)) & size(x,2)>0, Beginings,Ends);

    % length = end - start
    Dummy = cellfun(@(x,y) y-x, Beginings(GoodOnes),Ends(GoodOnes),'UniformOutput', false); 
    BurstLength(GoodOnes) = Dummy;
    
    F = mat2cell(MaxF,[size(MaxF,1)],[repmat(1,1,size(MaxF,2))]);
    

    % max F value in each burst
      %Dummy = cellfun(@(x,y,z) z(floor((y+x)./2))', Beginings(GoodOnes),Ends(GoodOnes),F(GoodOnes),'UniformOutput', false);
        Dummy = cellfun(@(x,y,z) arrayfun(@(a,b) nanmean(z(a:b)),x,y,'UniformOutput', false) , Beginings(GoodOnes),Ends(GoodOnes),F(GoodOnes),'UniformOutput', false);
    GoodIndex = find(GoodOnes);
    for i = 1:length(GoodIndex)
        BurstMax{GoodIndex(i)} = Cell2Mat(Dummy{i});
    end

    % position of max F value in each burst
    %Dummy = cellfun(@(x,y) (y+x)./2, Beginings(GoodOnes),Ends(GoodOnes),'UniformOutput', false);
        Dummy = cellfun(@(x,y,z) arrayfun(@(a,b) find(z(a:b) == max(z(a:b)),1)+a+StartCount+min(StartCount-1,1)-1,x,y,'UniformOutput', false), Beginings(GoodOnes),Ends(GoodOnes),F(GoodOnes),'UniformOutput', false);
    %BurstPeak(GoodOnes) = Dummy(:);
    for i = 1:length(GoodIndex)
        BurstPeak{GoodIndex(i)} = Cell2Mat(Dummy{i});
    end
    
     % burstsize, sum all NormF in each burst
    ElT = 5000/2000*60;
    Dummy = cellfun(@(x,y,z) arrayfun(@(a,b) nansum(z(a:b))* TimeRes ./ 2 ./ ElT,x,y,'UniformOutput', false) , Beginings(GoodOnes),Ends(GoodOnes),F(GoodOnes),'UniformOutput', false);
    GoodIndex = find(GoodOnes);
    for i = 1:length(GoodIndex)
        BurstSize{GoodIndex(i)} = Cell2Mat(Dummy{i});
    end
 
    
    % period = time from max to max
    Dummy = cellfun(@(x,y) diff([(x+y)./2]), Beginings(GoodOnes),Ends(GoodOnes),'UniformOutput', false);
    %Dummy = cellfun(@(x,y,z) diff(arrayfun(@(x,y) find(max(z(x:y))==z(x:y))+x,x,y,'ErrorHandler',@errorFunc), Beginings(GoodOnes),Ends(GoodOnes),F(GoodOnes),'UniformOutput', false);
    BurstPeriod(GoodOnes) = Dummy;
    
    % off time
    % only in traces that have more than 1 burst
    GoodOnes = cellfun(@(x,y) size(x,2)>1, Beginings,Ends);
    % off = time from on to on - the time that is actually on (off-on)
    %Dummy = cellfun(@(x,y) diff(x) - (y(1:end-1)-x(1:end-1)), Beginings(GoodOnes),Ends(GoodOnes),'UniformOutput', false);
    % time from one on (starting from 2nd) to the previous off (except last)
    Dummy = cellfun(@(x,y) x(2:end) - y(1:end-1), Beginings(GoodOnes),Ends(GoodOnes),'UniformOutput', false);
    OffTime(GoodOnes) = Dummy;
end