function[Data3D] = NuclearMembraneProperties(manual, Path, File, Name, Stats_tracked_3D, TimeRes, Frames, nc14, Delay, Data3D,minNumb)
%% MEASURE PROPERTIES in 2D and 3D
if manual
% get metadata file and select experiment to open
info = readtable(['~/Google Drive jf565/MATLAB_R_scripts/metadata MS2 3D',' Dillan','.txt'],'ReadVariableNames', true,'Delimiter', '\t');

Exp = 11; % experiment number to analyze
Table2Vars(info(Exp,:));

load([Path,File,Name,File,'_Stats3D.mat']);

end
PathToSave = [Path,File,Name,File]; 
load([Path,File,Name,File,'_Stats2D_Q1.mat']);
load([Path,File,Name,File,'_Stats2D_Q2.mat']);
load([Path,File,Name,File,'_Stats2D_Q3.mat']);
%Name = '_3D_nucmemb no imadjust, imdilate/'
%% extract properties from Stats and save in tables and Data3D structure (will save all stored properties)
%minNumb = 5
mkdir([Path,File,Name,'tables/'])

%[Data3D] = Stats2Struct(Stats_tracked_3D, TimeRes,Frames,minNumb, Path, File, Name,VarSuffix,FileSuffix);
[Data3D] = Stats2Struct(Data3D,Stats_tracked_3D, TimeRes,Frames,minNumb, Path, File, Name,'','3D');
[Data3D] = Stats2Struct(Data3D,Stats_tracked_2D_Q1, TimeRes,Frames,minNumb, Path, File, Name,'1','2DQ1');
[Data3D] = Stats2Struct(Data3D,Stats_tracked_2D_Q2, TimeRes,Frames,minNumb, Path, File, Name,'2','2DQ2');
[Data3D] = Stats2Struct(Data3D,Stats_tracked_2D_Q3, TimeRes,Frames,minNumb, Path, File, Name,'3','2DQ3');

 close all
        %%
        TimeScale = ([1:Frames] - nc14 + Delay).*TimeRes./60;
        Fig = figure('PaperOrientation','landscape');
subplot(451); plot(TimeScale,Data3D.Area1,'Color',[0.7,0.7,0.7]); hold on; plot(TimeScale,nanmean(Data3D.Area1,2),'.r'); title('Area')
subplot(452); plot(TimeScale,Data3D.Eccentricity1,'Color',[0.7,0.7,0.7]); hold on; plot(TimeScale,nanmean(Data3D.Eccentricity1,2),'.r'); title('Eccentricity')
subplot(453); plot(TimeScale,Data3D.Solidity1,'Color',[0.7,0.7,0.7]); hold on; plot(TimeScale,nanmean(Data3D.Solidity1,2),'.r'); title('Solidity')
subplot(454); plot(TimeScale,Data3D.Extent1,'Color',[0.7,0.7,0.7]); hold on; plot(TimeScale,nanmean(Data3D.Extent1,2),'.r');  title('Extent')
subplot(455); plot(TimeScale,Data3D.Perimeter1,'Color',[0.7,0.7,0.7]); hold on; plot(TimeScale,nanmean(Data3D.Perimeter1,2),'.r'); title('Perimeter')

subplot(456); plot(TimeScale,Data3D.Area2,'Color',[0.7,0.7,0.7]); hold on; plot(TimeScale,nanmean(Data3D.Area2,2),'.r')
subplot(457); plot(TimeScale,Data3D.Eccentricity2,'Color',[0.7,0.7,0.7]); hold on; plot(TimeScale,nanmean(Data3D.Eccentricity2,2),'.r')
subplot(458); plot(TimeScale,Data3D.Solidity2,'Color',[0.7,0.7,0.7]); hold on; plot(TimeScale,nanmean(Data3D.Solidity2,2),'.r')
subplot(459); plot(TimeScale,Data3D.Extent2,'Color',[0.7,0.7,0.7]); hold on; plot(TimeScale,nanmean(Data3D.Extent2,2),'.r')
subplot(4,5,10); plot(TimeScale,Data3D.Perimeter2,'Color',[0.7,0.7,0.7]); hold on; plot(TimeScale,nanmean(Data3D.Perimeter2,2),'.r')

subplot(4,5,11); plot(TimeScale,Data3D.Area3,'Color',[0.7,0.7,0.7]); hold on; plot(TimeScale,nanmean(Data3D.Area3,2),'.r')
subplot(4,5,12); plot(TimeScale,Data3D.Eccentricity3,'Color',[0.7,0.7,0.7]); hold on; plot(TimeScale,nanmean(Data3D.Eccentricity3,2),'.r')
subplot(4,5,13); plot(TimeScale,Data3D.Solidity3,'Color',[0.7,0.7,0.7]); hold on; plot(TimeScale,nanmean(Data3D.Solidity3,2),'.r')
subplot(4,5,14); plot(TimeScale,Data3D.Extent3,'Color',[0.7,0.7,0.7]); hold on; plot(TimeScale,nanmean(Data3D.Extent3,2),'.r')
subplot(4,5,15); plot(TimeScale,Data3D.Perimeter3,'Color',[0.7,0.7,0.7]); hold on; plot(TimeScale,nanmean(Data3D.Perimeter3,2),'.r')

subplot(4,5,16); plot(TimeScale,nanmean(Data3D.Area1,2)); hold on; plot(TimeScale,nanmean(Data3D.Area2,2));plot(TimeScale,nanmean(Data3D.Area3,2));
subplot(4,5,17); plot(TimeScale,nanmean(Data3D.Eccentricity1,2)); hold on; plot(TimeScale,nanmean(Data3D.Eccentricity2,2));plot(TimeScale,nanmean(Data3D.Eccentricity3,2));
subplot(4,5,18); plot(TimeScale,nanmean(Data3D.Solidity1,2)); hold on; plot(TimeScale,nanmean(Data3D.Solidity2,2));plot(TimeScale,nanmean(Data3D.Solidity3,2));
subplot(4,5,19); plot(TimeScale,nanmean(Data3D.Extent1,2)); hold on; plot(TimeScale,nanmean(Data3D.Extent2,2));plot(TimeScale,nanmean(Data3D.Extent3,2));
subplot(4,5,20); plot(TimeScale,nanmean(Data3D.Perimeter1,2)); hold on; plot(TimeScale,nanmean(Data3D.Perimeter2,2));plot(TimeScale,nanmean(Data3D.Perimeter3,2));

print(Fig,[PathToSave,'_2Dproperties.pdf'],'-fillpage', '-dpdf');

%%
Fig =  figure('PaperOrientation','landscape');
subplot(231); plot(TimeScale,Data3D.Volume,'Color',[0.7,0.7,0.7]); hold on; plot(TimeScale,nanmean(Data3D.Volume,2),'.r'); title('Volume')
subplot(232); plot(TimeScale,Data3D.Solidity,'Color',[0.7,0.7,0.7]); hold on; plot(TimeScale,nanmean(Data3D.Solidity,2),'.r'); title('Solidity')
subplot(233); plot(TimeScale,Data3D.SurfaceArea,'Color',[0.7,0.7,0.7]); hold on; plot(TimeScale,nanmean(Data3D.SurfaceArea,2),'.r'); title('SurfaceArea')
subplot(234); plot(TimeScale,Data3D.PrincipalAxisLength_1,'Color',[0.7,0.7,0.7]); hold on; plot(TimeScale,nanmean(Data3D.PrincipalAxisLength_1,2),'.r'); title('Principal Axis Length 1')
subplot(235); plot(TimeScale,Data3D.PrincipalAxisLength_2,'Color',[0.7,0.7,0.7]); hold on; plot(TimeScale,nanmean(Data3D.PrincipalAxisLength_2,2),'.r'); title('Principal Axis Length 2')
subplot(236); plot(TimeScale,Data3D.PrincipalAxisLength_3,'Color',[0.7,0.7,0.7]); hold on; plot(TimeScale,nanmean(Data3D.PrincipalAxisLength_3,2),'.r'); title('Principal Axis Length 3')

print(Fig,[PathToSave,'_3Dproperties.pdf'],'-fillpage', '-dpdf');
end