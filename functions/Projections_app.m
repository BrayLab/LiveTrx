 function[]=Projections_app(Parameters)

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

RXYT = zeros(Height, Width,Frames);
GXYT = zeros(Height, Width,Frames);
RXZT = zeros(Slices,Width,Frames);
GXZT = zeros(Slices,Width,Frames);

if Channels == 3; 
    TLend = zeros(Height, Width,Frames);
end
%
for f = 1:Frames
    disp(['reading f',num2str(f),'...']);
    G = ReadSingleStack(reader,Channels,Slices,Frames0,Flip,From,1,f);
    GXYT(:,:,f) = max(G,[],3);
    GXZT(:,:,f) = permute(max(G,[],1),[3,1,2]);
    
    if Channels >1
        R = ReadSingleStack(reader,Channels,Slices,Frames0,Flip,From,2,f);
        RXYT(:,:,f) = max(R,[],3);
        RXZT(:,:,f) = permute(max(R,[],1),[3,1,2]);
    end
    
    if Channels == 3
       TL = ReadSingleStack(reader,Channels,Slices,Frames0,Flip,From,3,f);
       if strcmp(Channel3,'TL')
            TLend(:,:,f) = TL(:,:,Slices);
       else
            TLend(:,:,f) = max(TL,[],3);   
       end
           
    end
end
%
if Channels > 1
    [XYT] = CombineGR(GXYT,RXYT,Bits);
    [XZT] = CombineGR(GXZT,RXZT,Bits);
    Combined = cat(1,XZT,XYT);
else
    [XYT] = CombineGR(GXYT,[],Bits);
    [XZT] = CombineGR(GXZT,[],Bits);
    Combined = cat(1,XZT,XYT);
end
% time stamp
%if ~manual
    [Combined] = TimeStamp(Combined,TimeRes,nc14,Delay);
%end
WriteRGB(Combined, PathToSave, '_CombinedY.tiff','none')
%
if Channels == 3
    TLRGB = permute(cat(4,TLend,TLend,TLend),[1,2,4,3]).*255./(2^Bits-1);
    Combined = cat(2,XYT,TLRGB);
    [Combined] = TimeStamp(Combined,TimeRes,nc14,Delay);
    WriteRGB(Combined, PathToSave, '_CombinedTL.tiff','none')
end
% time stamp

    GXY = max(GXYT,[],3);
    GXT = permute(max(GXYT,[],1),[3,2,1]);
    GYT = permute(max(GXYT,[],2),[1,3,2]);
    GZT = permute(max(GXZT,[],2),[1,3,2]);
if Channels > 1
    RXY = max(RXYT,[],3);
    RXT = permute(max(RXYT,[],1),[3,2,1]);
    RYT = permute(max(RXYT,[],2),[1,3,2]);
    RZT = permute(max(RXZT,[],2),[1,3,2]);
    
    [XY] = CombineGR(GXY,RXY,Bits);
    [XT] = CombineGR(GXT,RXT,Bits);
    [YT] = CombineGR(GYT,RYT,Bits);
    [ZT] = CombineGR(GZT,RZT,Bits);
else
    [XY] = CombineGR(GXY,[],Bits);
    [XT] = CombineGR(GXT,[],Bits);
    [YT] = CombineGR(GYT,[],Bits);
    [ZT] = CombineGR(GZT,[],Bits);
end
 imwrite(XY,jet, [PathToSave, '_XY.tiff'],'Compression','none')
imwrite(XT,jet, [PathToSave, '_XT.tiff'],'Compression','none')
imwrite(YT,jet, [PathToSave, '_YT.tiff'],'Compression','none')
imwrite(ZT,jet, [PathToSave, '_ZT.tiff'],'Compression','none')
