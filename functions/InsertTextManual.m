function[Im] = InsertTextManual(Im,PosXY,Text,FontSize,Grey)
    %PosXY = [50,10]; Text ='test'; FontSize = 14;Grey = 255;
    
    Fig = figure('color','black','visible','off','position',[1 1 size(Im,2) size(Im,1)]); 
    imshow(zeros(size(Im))); 
    set(gca,'units','pixels','position',[1 1 size(Im,2) size(Im,1)])
    text('units','pixels','position',PosXY,'fontsize',FontSize,'string',Text,'Color','white') 
    Frame = getframe; 
    close(Fig) 
    %Frame = rgb2gray(Frame.cdata);
    %Mask = (Frame ~= 0); 
    %Im(repmat(Mask,1,1,3)) = uint8(Grey); 
    Im = Im.*uint8(~Frame.cdata) + Frame.cdata;
   
end