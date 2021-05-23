function[Bits,Width,Height, Channels, Slices, Frames,Frames0, XRes, YRes, ZRes,Zoom,TimeRes, Notes] = readMetadataBFOME(reader)


%%
%[Series, OME, Org] = GetMetadata(reader)
%globalMetadata = string(reader.getCoreMetadata()); %commented 18/07/20
%InfoGlobal = string(transpose(strsplit(globalMetadata{1},'\n'))); %commented 18/07/20
%globalMetadata = reader.getGlobalMetadata(); %empty
%reader.setSeries(2);
%reader.getSeriesMetadataValue();
seriesMetadata = reader.getSeriesMetadata(); 
%javaMethod('merge', 'loci.formats.MetadataTools', ...
%               globalMetadata, seriesMetadata, 'Global ');
%           javaMethod('getGlobalMetadata', reader)
   %           javaMethod('getCoreMetadataList', reader)
%javaMethod('getMetadataValue', reader)
%Info = vertcat(InfoGlobal, strsplit(string(seriesMetadata),', ')'); %commented 18/07/20
Info = strsplit(string(seriesMetadata),', ')'; %uncommented 18/07/20
Bits = FindInMetadata(Info,'BitSize');
TimeRes = FindInMetadata(Info,'CycleTime');
Pinhole = FindInMetadata(Info,'|PinholeAiry');
ScanMode = FindInMetadata(Info,'ScanMode');
ScanSpeed = FindInMetadata(Info,'ScanSpeed');
Zoom = FindInMetadata(Info,'|Zoom');
Order = FindInMetadata(Info,'DimensionOrder');
%Width = FindInMetadata(Info,'SizeX');
%Height = FindInMetadata(Info,'SizeY');
%Slices = FindInMetadata(Info,'SizeZ');
%Channels = FindInMetadata(Info,'SizeC');
%Frames = FindInMetadata(Info,'SizeT');
GreenP = FindInMetadata(Info,'|IntensityDev 3');
RedP = FindInMetadata(Info,'|IntensityDev 6');
GreenG = FindInMetadata(Info,'|Gain 1');
RedG = FindInMetadata(Info,'|Gain 3');
A488P = FindInMetadata(Info,'|OutputPowerPercentage');
A488W = FindInMetadata(Info,'|OutputPowerWatt 1');
Begin = FindInMetadata(Info,'|Begin');
End = FindInMetadata(Info,'|End');

Notes = ['Pinh_', num2str(round(Pinhole,2)),'_Speed_',num2str(ScanSpeed),'_488_',num2str(round(A488P,2)),'_G_',num2str(round(GreenP,1)), '_R_', num2str(round(RedP,2)),'_GG_',num2str(round(GreenG,1)), '_RG_', num2str(round(RedG,1))];
% Index = find(cellfun(@(x) isempty(x), strfind(Info,'BitSize')) == 0);
% D = Info{Index}; D = strsplit(D,'='); Bits = str2num(D{2});
% Index = find(cellfun(@(x) isempty(x), strfind(Info,'CycleTime')) == 0);
% D = Info{Index}; D = strsplit(D,'='); TimeRes = str2num(D{2});
% Index = find(cellfun(@(x) isempty(x), strfind(Info,'|PinholeAiry')) == 0);
% D = Info{Index}; D = strsplit(D,'='); Pinhole = str2num(D{2});
% Index = find(cellfun(@(x) isempty(x), strfind(Info,'ScanMode')) == 0);
% D = Info{Index}; D = strsplit(D,'='); ScanMode = D{2};
% Index = find(cellfun(@(x) isempty(x), strfind(Info,'ScanSpeed')) == 0);
% D = Info{Index}; D = strsplit(D,'='); ScanSpeed = str2num(D{2});
% Index = find(cellfun(@(x) isempty(x), strfind(Info,'|Zoom')) == 0);
% D = Info{Index}; D = strsplit(D,'='); Zoom = str2num(D{2});
%%
omeMeta = reader.getMetadataStore();
Info = strsplit(char(omeMeta.dumpXML()),' ')';
XRes = extractInfo(Info, '^PhysicalSizeX=');XRes = mean(XRes);
YRes = extractInfo(Info, '^PhysicalSizeY=');YRes = mean(YRes);
ZRes = extractInfo(Info, '^PhysicalSizeZ=');ZRes = mean(ZRes);
Width = extractInfo(Info, '^SizeX='); Width = mean(Width);
Height = extractInfo(Info, '^SizeY='); Height = mean(Height);
Slices = extractInfo(Info, '^SizeZ='); Slices = mean(Slices); 
Channels = extractInfo(Info, '^SizeC='); Channels = mean(Channels);
Frames0 = extractInfo(Info, '^SizeT=');
Frames = sum(Frames0);
if isnan(ZRes)==1
    ZRes = (End - Begin)./(Slices-1)*10^6;
end
end
