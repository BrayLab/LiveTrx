function[] = PlotsPosition(OnOff,TimeScale,Selected,Properties,CentX,CentY,CentZ,SpotZ,XYRes,ZRes,Width, Height, FileOut)
Fig = figure('PaperSize',[50 50],'PaperUnits','inches','resize','on', 'visible','off');
Fig.Renderer='Painters';
OnOffSelected = OnOff(:,Selected);
PropertiesSelected = Properties(Selected,:);
CentXSelected = CentX(:,Selected).*XYRes;
CentYSelected = CentY(:,Selected).*XYRes;
CentZSelected = CentZ(:,Selected).*ZRes;
subplot(221);plot(TimeScale,CentXSelected); hold on; xlabel('t(min)'); ylabel('Xpos (um)')
plot(TimeScale,nanmean(CentXSelected,2),'.r','MarkerSize',12);
subplot(222);plot(TimeScale,CentYSelected); hold on; xlabel('t(min)'); ylabel('Ypos (um)')
plot(TimeScale,nanmean(CentYSelected,2),'.r','MarkerSize',12);
%subplot(223);plot(TimeScale,CentZSelected); hold on; xlabel('t(min)'); ylabel('Zpos (um)')
%plot(TimeScale,nanmean(CentZSelected,2),'.r','MarkerSize',12);
subplot(223);plot(CentXSelected,CentYSelected); xlabel('X (um))'); ylabel('Y (um)');
xlim([1, Width].*XYRes); ylim([1,Height].*XYRes);
set (gca,'Ydir','reverse')

CentZSelected = CentZ(:,Selected).*OnOffSelected.*ZRes;CentZSelected(CentZSelected==0)= NaN;
SpotZSelected = SpotZ(:,Selected).*OnOffSelected.*ZRes;SpotZSelected(SpotZSelected==0)= NaN;
DistZSelected = SpotZSelected-CentZSelected;
subplot(224);
errorbar(TimeScale,nanmean(CentZSelected,2),nanstd(CentZSelected,0,2)./sum(~isnan(CentZSelected),2),'.r'); hold on
errorbar(TimeScale,nanmean(SpotZSelected,2),nanstd(CentZSelected,0,2)./sum(~isnan(SpotZSelected),2),'.g'); 
errorbar(TimeScale,nanmean(DistZSelected,2),nanstd(DistZSelected,0,2)./sum(~isnan(DistZSelected),2),'.b');
plot(TimeScale, TimeScale.*0);
xlabel('t(min)'); ylabel('Zpos (um)')
set (gca,'Ydir','reverse')
print(Fig,FileOut,'-fillpage', '-dpsc','-append');

end