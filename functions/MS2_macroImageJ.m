 function[]=MS2_macroImageJ(PathToSave,manual)
% BEFORE START DO manual = true
%%
if manual == true
    [File,Path] = uigetfile('*.*','select file for name / MAX.tif','/Volumes/Mac OS/MS2/');
    Name = '_3D/'
    Channel1 = 'MCP'; Channel2 = 'his';
    Nickname='m5m8peve mDr:TTG'; Rep = 4;
        Nickname='m5m8peveIII kukRi m&z'; Rep = 3
                 Nickname='m5m8peve Fog'; Rep = 4
                         Nickname='m5m8peve kukPE'; Rep = 3; SplitEarly = 15;
                         %Nickname='m5m8peveIII Nullo'; Rep = 3; SplitEarly = 15;

                %Nickname='m8NEpevex2'; Rep = 3
                %Nickname='m5m8Dtwidlpevex2'; Rep = 5

        %Nickname='vndEEEpeve ctaRi'; Rep = 7
                %Nickname='vndEEEpeve Fog:+'; Rep = 5; SplitEarly = 0;
        %Nickname='m5m8peve + eveNICD'; Rep = 5; SplitEarly = 15;
        %Nickname='pUbi'; Rep = 1; SplitEarly = 0;
        %Nickname='nup107RFP'; Rep = 1  
            
            MetaFile = ' enhprom';
            %MetaFile = ' ecNICD';
            MetaFile = ' mutBG';
            %MetaFile = ' other';
            %MetaFile = ' javier';
           % MetaFile = ' 2c';
           %MetaFile = '';
    PathToSave = [Path,File];
end  
    %PathToSave = [Path,File]; 
    Flip = [0,0] % first to flip vertically 2nd to flip horizontally
    From = 1;%230
    nc14 = 1;
    To='NA'; %NA
    Delay=0;
    Notes = '';
    try
        [Dummy,Comments] = system(['mdls -raw -name kMDItemFinderComment "',Path,File,'"'])
        Comments = strsplit(Comments,', ');
        try;Index = find(cellfun(@(x) isempty(x), strfind(lower(Comments),'flip')) == 0);
        D = Comments{Index}; D = strsplit(D,' '); Flip = [str2num(D{2}(1)),str2num(D{2}(2))];end
        try;Index = find(cellfun(@(x) isempty(x), strfind(lower(Comments),'nc14')) == 0);
        D = Comments{Index}; D = strsplit(D,' '); nc14 = str2num(D{2});end
        try;Index = find(cellfun(@(x) isempty(x), strfind(lower(Comments),'delay')) == 0);
        D = Comments{Index}; D = strsplit(D,' '); Delay = str2num(D{2});end
        try;Index = find(cellfun(@(x) isempty(x), strfind(lower(Comments),'to')) == 0);
        D = Comments{Index}; D = strsplit(D,' '); To = str2num(D{2});end
    catch
        disp('Comments not found! set Flip, nc14 and delay');   
    end

reader = bfGetReader([PathToSave])
[Bits,Width,Height, Channels, Slices, Frames, Frames0, XRes, YRes, ZRes,Zoom,TimeRes,Settings] = readMetadataBFOME(reader);
try
    XYRes=mean(round(mean([XRes,YRes]),2));
    ZRes = mean(round(ZRes,2));
end

if strcmp(To,'NA')==1 | To == inf;  To=Frames-1;
else
    if isnumeric(To) == 0 
        To = str2num(To); end; end;
Frames = To-From+1;
%
%

if manual

Metadata = readtable(['~/Google Drive jf565/MATLAB_R_scripts/metadata MS2 3D',MetaFile,'.txt'],'Delimiter', '\t');

NewMetadata = cell2table({Path,File,Name,Frames,Slices,Bits,TimeRes,XYRes,ZRes,...
  Delay,Nickname, Rep,Zoom,strjoin(arrayfun(@(x) num2str(x),Flip,'UniformOutput',false),','),From,nc14,To,SplitEarly,Channel1,Channel2,[Notes,' ',Settings]},...
    'VariableNames', {'Path','File','Name','Frames','Slices', 'Bits','TimeRes','XYRes','ZRes','Delay','Nickname','Rep','Zoom','Flip','From','nc14','To','SplitEarly','Channel1','Channel2','Notes'});
SaveMetadata = [Metadata;NewMetadata];
writetable(SaveMetadata,['~/Google Drive jf565/MATLAB_R_scripts/metadata MS2 3D',MetaFile,'.txt'],'Delimiter', '\t');
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
       TLend(:,:,f) = TL(:,:,Slices);
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
if ~manual
    [Combined] = TimeStamp(Combined,TimeRes,nc14,Delay);
end
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
