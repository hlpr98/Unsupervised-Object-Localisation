function WindowButtonMotionCallback(a, b)

global h I pcs proposals selections current

%% Ask for center:
p = get(h, 'currentPoint');
H = size(I, 1);
W = size(I, 2);
i = round(H-p(2));
j = round(p(1));
id = 1;

%% Find closest proposal:
d = sqrt((pcs(:, 1) - j) .^ 2 + (pcs(:, 2) - i) .^ 2);
[~, ids] = min(d);
current = proposals(ids(id), :);

%% Display results:
%clf;

% windowSize = 5;
% % kernel = ones(windowSize) / windowSize^2; % Mean Blur
% % kernel = [1,2,1;2,4,2;1,2,1]/windowSize^2;  % Gaussian Blur -- windlowSize = 3x3
% kernel = [1,1,1,1,1;1,3,5,3,1;1,5,10,5,1;1,3,5,3,1;1,1,1,1,1]/windowSize^2;  % Gaussian Blur -- windlowSize = 5x5
% 
% grayImage = rgb2gray(I);
% [Gmag,Gdir] = imgradient(grayImage);
% blurIm = conv2(double(Gmag), kernel, 'same');
% blurIm = uint8(blurIm);
% [R,C] = size(blurIm);
% for n = 1:R
%     for m = 1:C
%         if blurIm(n,m)<250
%             blurIm(n,m) = 0;
%         end
%     end
% end    
% % edgeIm = edge(blurIm, 'canny');
% imshow(blurIm, 'Border', 'tight');


imshow(I, 'Border', 'tight');
hold on;
plotBoxes(selections, 'b', [], '-');
plotBoxes(current, 'y', [], '-');
plot([pcs(ids(id), 1) j], [pcs(ids(id), 2) i], 'g-.');
plot(pcs(ids(id), 1),pcs(ids(id), 2), 'yx', 'LineWidth', 5);
plot(j, i, 'bx', 'LineWidth', 5);

end