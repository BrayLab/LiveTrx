function[Stack3D] = ReadSingleStack(reader,Channels,Slices,Frames0,Flip,From,Channel,f)
    CFrames = cumsum(Frames0);
    Series = find(f+From <= CFrames & f+From > CFrames-Frames0 == 1);
    reader.setSeries(Series - 1);
    Frame= f+From-CFrames(Series)+Frames0(Series); c = Channel; z = 1:Slices;
    k = (Frame-1)*Channels*Slices + (z-1)*Channels + c;
    SingleFrame = arrayfun(@(x) double(bfGetPlane(reader,x)),k, 'UniformOutput', false);
    
    for i = 1:size(SingleFrame,2)
        Stack3D(:,:,i) = SingleFrame{i};
    end
%     if Flip(1) == 1
%     	Stack3D = flip(Stack3D,1);
%     end
%     if Flip(2) == 1
%         Stack3D = flip(Stack3D,2);
%     end
%     try Flip(3);
%             if Flip(3) == 1
%                 Stack3D = flip(Stack3D,3);
%             end
%     end
    IndexFlip = find(Flip);
    for i = 1:length(IndexFlip)
        Stack3D = flip(Stack3D,IndexFlip(i));
    end
end