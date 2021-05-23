function[RGB] = CombineGR(Green,Red,Bits)
if isempty(Red)
    Red = zeros(size(Green));
end
for f = 1:size(Red,3)
    RGB1(:,:,1,f) = uint8(Green(:,:,f).*255./(2^Bits-1));
    RGB1(:,:,2,f) = uint8(Green(:,:,f).*255./(2^Bits-1));
    RGB1(:,:,3,f) = uint8(Green(:,:,f).*255./(2^Bits-1));
    RGB2(:,:,1,f) = uint8(zeros(size(Red,1),size(Red,2)));
    RGB2(:,:,2,f) = uint8(Green(:,:,f).*255./(2^Bits-1));
    RGB2(:,:,3,f) = uint8(Red(:,:,f).*255./(2^Bits-1));
    RGB(:,:,:,f) = cat(2, RGB1(:,:,:,f), RGB2(:,:,:,f));
end
end