function[] = Write16b(toWrite, PathToSave, Suffix)
imwrite(uint8(toWrite(:,:,1)), [PathToSave, Suffix],'Compression','none')
for f = 2:size(toWrite,3)
    imwrite(uint8(toWrite(:,:,f)), [PathToSave, Suffix],'WriteMode','append','Compression','none')
end
end