function[Combined] = TimeStamp(Combined,TimeRes,nc14,Delay)
    TimeScale = TimeRes*(-nc14+Delay-1);
    %Fig = figure('color','white','visible','off'); 
    for f = 1:size(Combined,4)
        disp(['f',num2str(f)])
        TimeScale = TimeScale + TimeRes;
        if TimeScale < 0; Sign = '-'; else; Sign =' ';end
        Hour = abs(fix(TimeScale./3600));
        Min = fix(abs(rem(TimeScale,3600))./60);
        Sec = rem(abs(rem(TimeScale,3600)),60);
        FormatSpec = '%0d:%02d:%02d';
        Text = [Sign,sprintf(FormatSpec,Hour,Min,Sec)];
        %Text = datestr(datenum(minutes(floor(TimeScale./60))+seconds(mod(TimeScale,60))),'Format','mm:ss')
        try
            Combined(:,:,:,f) = insertText(Combined(:,:,:,f),[size(Combined,2)-100,size(Combined,1)-30],Text,'FontSize',20,'TextColor','white','BoxOpacity',0);
        catch
            Combined(:,:,:,f) = InsertTextManual(Combined(:,:,:,f),[size(Combined,2)-100,30],Text,20,255);
            %[Combined(:,:,:,f),Fig] = InsertTextManual2(Fig,Combined(:,:,:,f),[size(Combined,2)-100,30],Text,20,255);
        end
    end
    %close(Fig)
end