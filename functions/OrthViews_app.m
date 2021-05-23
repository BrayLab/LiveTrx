 function[]=OrthViews_app(Parameters)

 Table2Vars(Parameters);
Flip = str2double(strsplit(Flip,','));

mkdir([Path,File,Name]) %uncomment
PathToSave = [Path,File]; %uncomment
show = 'off';
reader = bfGetReader([Path,File]);
[Bits,Width,Height, Channels, Slices, Frames, Frames0, XRes, YRes, ZRes,Zoom,TimeRes,Settings] = readMetadataBFOME(reader);
try
    XYRes=mean(round(mean([XRes,YRes]),2));
    ZRes = mean(round(ZRes,2));
end
if strcmp(To,'NA')==1 | To == inf;  To=Frames-1;
elseif isnumeric(To) == 0; To = str2num(To); end;
Frames = To-From+1;

if strcmp(Channel1,'nucmemb') | strcmp(Channel2,'nucmemb')
    NucMemb = 1;
else
    NucMemb = 0
end
if strcmp(Channel1,'memb') | strcmp(Channel2,'memb')
    Membranes = 1;
else
    Membranes = 0
end


%%
NewSlices = round(Slices*ZRes/XYRes);
RXYT = zeros(Height,Width,Frames);
GXYT = zeros(Height,Width,Frames);
RXZT = zeros(NewSlices,Width,Frames);
GXZT = zeros(NewSlices,Width,Frames);
RYZT = zeros(Height,NewSlices,Frames);
GYZT = zeros(Height,NewSlices,Frames);

if Channels == 3; 
    TLend = zeros(Height, Width,Frames);
end
%
StartSliceMembProj = 5;
StartSliceMembProj = 10;
OrthDepth = 3;
for f = 1:Frames
    disp(['reading f',num2str(f),'...']);
    
    if strcmp(Channel1,'memb') | strcmp(Channel1,'his') | strcmp(Channel1,'nucmemb')
        G = ReadSingleStack(reader,Channels,Slices,Frames0,Flip,From,1,f);
        GXYT(:,:,f) = max(G(:,:,StartSliceMembProj:Slices-5),[],3);
        GYZT(:,:,f) = imresize(permute(max(G(:,round(Width/2)-OrthDepth:round(Width/2)+OrthDepth,:),[],2),[1,3,2]),[Height,NewSlices]);
        %GYZT(:,:,f) = imresize(permute(max(G(:,:,:),[],2),[1,3,2]),[Height,NewSlices]);
        GXZT(:,:,f) = imresize(permute(max(G(round(Height/2)-OrthDepth:round(Height/2)+OrthDepth,:,:),[],1),[3,2,1]),[NewSlices,Width]);
    end
    if strcmp(Channel1,'MCP') 
       G = ReadSingleStack(reader,Channels,Slices,Frames0,Flip,From,1,f);
        GXYT(:,:,f) = max(G,[],3);
        GYZT(:,:,f) = imresize(permute(max(G,[],2),[1,3,2]),[Height,NewSlices]);
        GXZT(:,:,f) = imresize(permute(max(G,[],1),[3,2,1]),[NewSlices,Width]);
    end
    
    if Channels > 1 & (strcmp(Channel2,'memb') | strcmp(Channel2,'his')| strcmp(Channel2,'nucmemb') )      
        R = ReadSingleStack(reader,Channels,Slices,Frames0,Flip,From,2,f);
        RXYT(:,:,f) = max(R(:,:,StartSliceMembProj:Slices-5),[],3);
        RYZT(:,:,f) = imresize(permute(max(R(:,round(Width/2)-OrthDepth:round(Width/2)+OrthDepth,:),[],2),[1,3,2]),[Height,NewSlices]);
        %RYZT(:,:,f) = imresize(permute(max(R(:,:,:),[],2),[1,3,2]),[Height,NewSlices]);
        RXZT(:,:,f) = imresize(permute(max(R(round(Height/2)-OrthDepth:round(Height/2)+OrthDepth,:,:),[],1),[3,2,1]),[NewSlices,Width]);
    elseif Channels > 1
        R = ReadSingleStack(reader,Channels,Slices,Frames0,Flip,From,2,f);
        RXYT(:,:,f) = max(R,[],3);
        RYZT(:,:,f) = imresize(permute(max(R,[],2),[1,3,2]),[Height,NewSlices]);
        RXZT(:,:,f) = imresize(permute(max(R,[],1),[3,2,1]),[NewSlices,Width]);        
    end
    
%     if Channels == 3
%        TL = ReadSingleStack(reader,Channels,Slices,Frames0,Flip,From,3,f);
%        TLend(:,:,f) = TL(:,:,Slices);
%     end
end
%
% if Channels > 1
%     [XYT] = CombineGR(GXYT,RXYT,Bits);
%     [XZT] = CombineGR(GXZT,RXZT,Bits);
%     Combined = cat(1,XZT,XYT);
% else
%     [XYT] = CombineGR(GXYT,[],Bits);
%     [XZT] = CombineGR(GXZT,[],Bits);
%     Combined = cat(1,XZT,XYT);
% end
%% measure membrane length

% if strcmp(Channel1,'memb')
%     try
%         MaxHeight1 = readtable([Path,File,'_MembraneLengthCh1.txt'])
%     catch
%         [MaxHeight1] = MeasureMembranes(GXZT,Bits);
%         writetable(array2table(MaxHeight1'),[Path,File,'_MembraneLengthCh1.txt'],'Delimiter', '\t'); 
%     end
% end
% 
% if Channels > 1 & strcmp(Channel2,'memb')   
%     try
%         MaxHeight2 = readtable([Path,File,'_MembraneLengthCh2.txt'])
%     catch
%         [MaxHeight2] = MeasureMembranes(RXZT,Bits);
%         writetable(array2table(MaxHeight2'),[Path,File,'_MembraneLengthCh2.txt'],'Delimiter', '\t');  
%     end
% end

%% combine and save

if Channels >= 2
    Combined1 = cat(2,zeros(NewSlices,NewSlices,Frames),RXZT,GXZT);
    Combined2 = cat(2,RYZT,RXYT,GXYT);
    Combined = cat(1,Combined1,Combined2).*2^(8-Bits);
    Combined = cat(4,Combined,Combined,Combined);
    Combined = permute(Combined,[1,2,4,3]);
end

if Channels == 1
    Combined1 = cat(2,zeros(NewSlices,NewSlices,Frames),GXZT);
    Combined2 = cat(2,GYZT,GXYT);
    Combined = cat(1,Combined1,Combined2).*2^(8-Bits);
    Combined = cat(4,Combined,Combined,Combined);
    Combined = permute(Combined,[1,2,4,3]);
end

    [Combined] = TimeStamp(uint8(Combined),TimeRes,nc14,Delay);

WriteRGB(Combined, PathToSave, '_CombinedOrthViews.tiff','none')
%
% if Channels == 3
%     TLRGB = permute(cat(4,TLend,TLend,TLend),[1,2,4,3]).*255./(2^Bits-1);
%     Combined = cat(2,XYT,TLRGB);
%     [Combined] = TimeStamp(Combined,TimeRes,nc14,Delay);
%     WriteRGB(Combined, PathToSave, '_CombinedTL.tiff','none')
end

