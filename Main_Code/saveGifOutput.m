function saveGifOutput(F1, videoName)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% frame = F1(1);
% im = frame2im(frame);
% [imind, cm] = rgb2ind(im, 256);
% imwrite(imind, cm, videoName, 'gif', 'LoopCount', inf);

for i = 1:length(F1)
    frame = F1(i);
    im = frame2im(frame);
    [imind, cm] = rgb2ind(im, 256);
    imwrite(imind, cm, videoName, 'gif', 'WriteMode', 'append');
end
