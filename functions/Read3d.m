%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% lam94@cam.ac.uk, 2013
% Read tif stack with name filename
% Input: filename
% Output:   A   - 3d matrix with image data
%           info - information on the image
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [A info]=  Read3d(filename);
    info = imfinfo(filename,'tif');
    num_images = numel(info);   
    A = zeros(info(1).Height(1),info(1).Width ,num_images);
    for k = 1:num_images
            A(:,:,k) = double(imread(filename, k, 'Info', info));
    end
end