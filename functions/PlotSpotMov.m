function[] = PlotSpotMov(FileOut, Spots, Selected,Properties, CentX, CentY, SpotX, SpotY,varargin)
    CentX = CentX(:,Selected);
    CentY = CentY(:,Selected);
    SpotX = SpotX(:,Selected);
    SpotY = SpotY(:,Selected);
    PropertiesSelected = Properties(Selected,:);
    XPos = SpotX - CentX;
    YPos = SpotY - CentY;
    if length(varargin) > 0
         SpotXGauss = varargin{1};    SpotXGauss = SpotXGauss(:,Selected);
         SpotYGauss = varargin{2};    SpotYGauss = SpotYGauss(:,Selected);
         XPos1 = SpotX - SpotXGauss;
         YPos1 = SpotY - SpotYGauss;
        if length(varargin) > 2
            SpotXGauss2 = varargin{3};    SpotXGauss2 = SpotXGauss2(:,Selected);
            SpotYGauss2 = varargin{4};    SpotYGauss2 = SpotYGauss2(:,Selected);
            XPos2 = SpotX - SpotXGauss2;
            YPos2 = SpotY - SpotYGauss2;
        end
    end 
    
    
    for n = 1:(floor(length(find(Selected))/30)+1)
    disp([num2str(n),'/',num2str((floor(length(find(Selected))/30)+1))])
    Fig = figure('PaperSize',[30 50],'PaperUnits','inches','resize','on', 'visible','off');
    Ncol = 4; Nrow = 8;
        for i=1:Ncol*Nrow
            try
            subplot(Nrow,Ncol,i); 
            plot(XPos(:,(n-1)*Ncol*Nrow+i),YPos(:,(n-1)*Ncol*Nrow+i),'.-','Color',[0.7,0.7,0.7,0.2]); axis equal; hold on; % max pixel value
            if Spots > 0
                plot(XPos1(:,(n-1)*Ncol*Nrow+i),YPos1(:,(n-1)*Ncol*Nrow+i),'.-','Color',[54,100,139]./255); % gauss spot 1
                if Spots == 2
                    plot(XPos2(:,(n-1)*Ncol*Nrow+i),YPos2(:,(n-1)*Ncol*Nrow+i),'.-','Color',[105,139,34]./255); % gauss spot 2
                end
            end
            title(join(['#',num2str(PropertiesSelected.Label((n-1)*Ncol*Nrow+i)),' ',PropertiesSelected.Type((n-1)*Ncol*Nrow+i)]));
            hold off
            end
        end
     if n==1
          print(Fig,FileOut,'-fillpage', '-dpsc');
     else
          print(Fig,FileOut,'-fillpage', '-dpsc','-append');
     end
    end
    close all
end
    