% this function runs the PIV
% INPUT:    *	movie directory and .tif file [directory, file]
%           *	iteration [k]
%           *	PIV parameters [params]
%           *	variable stating if the PIV is run in full or just on the first
%               frame (test to correlate frame) [happy_test]
% OUTPUT:   *   raw PIV field [vraw]
%           *   .tif stack of the raw vector field (only if running full
%               PIV on all frames)

function [x_p, y_p, x_v, y_v, c_val] = run_piv(directory, file, k, params, happy_test, output_name)

%% IMAGE OUTLINE

% load current and next image
im1 = double(imread(fullfile(directory, file), k)) / 255; % transform image to double
im2 = double(imread(fullfile(directory, file), k + 1)) / 255;

% remove small particles and correct shape (erosion/dilation)
BW1 = edge(im1, 'Canny');	% find edges
dilate_radius = 4;
connectivity_fill = 4;
erode_radius = 4;
BW2d = imdilate(BW1, strel('disk', dilate_radius));
BW2f = imfill(BW2d, connectivity_fill, 'holes');
BW2 = imerode(BW2f, strel('disk', erode_radius));

[L,~] = bwlabel(BW2);

stats = regionprops(L, 'Area');
allArea = [stats.Area];
area_largest_obj = max(allArea(:));

BW2 = bwareaopen(BW2,area_largest_obj);

B = bwboundaries(BW2);
y = smooth(B{1,1}(:,1),6);
x = smooth(B{1,1}(:,2),6);

eroded_BW2 = imerode(BW2, strel('disk', round(params.search_size/2)));
% eroded_BW2 = imerode(BW2, strel('disk', erode_radius));
%% PIV

n = 1;

for i =  1:params.grid_distance:size(im1, 1)
    
    for j = 1:params.grid_distance:size(im1, 2)
        
        if eroded_BW2(i,j) == 1
            
            try
                sub_roi = im1(i-params.source_size:i+params.source_size,...
                    j-params.source_size:j+params.source_size);
            catch
                continue
            end
            
            try
                sub_area = im2(i-params.search_size:i+params.search_size,...
                    j-params.search_size:j+params.search_size);
            catch
                continue
            end
            
            try
                c = normxcorr2(sub_roi, sub_area);
            catch
                continue
            end
            
            [xroi, yroi] = size(sub_roi);
            [xarea, yarea] = size(sub_area);
            c = c(xroi:xarea, yroi:yarea);
            [max_c, imax] = max(c(:));
            [ypeak, xpeak] = ind2sub(size(c), imax(1));
            [xc, yc] = size(c);
            n = n + 1;
            
            if max_c >= params.correlation_threshold
                x_p(n) = j;
                y_p(n) = i;
                x_v(n) = xpeak - (xc + 1) / 2;
                y_v(n) = ypeak - (yc + 1) / 2;
                c_val(n) = max_c;
            end
            
        end
    end
end

%% VISUALISATION



% show current frame
imshow(im1, []);

% add cell edge in magenta
hold on
plot(x,y,'Color','m','LineWidth',1.5);
drawnow

% add vector field in green
quiver(x_p, y_p, x_v, y_v,'g','LineWidth',1);
drawnow
hold off
pause(1)

%% OUTPUT

% % save vector field in [vraw]
% pause(1)
% vraw(k).x = x_p;
% vraw(k).y = y_p;
% vraw(k).vx = x_v;
% vraw(k).vy = y_v;
% vraw(k).cc = c_val;
% clear x_p y_p x_v y_v c_val

% save figures in [images] if running the PIV on all frames
if happy_test == 1
    im_vectors = getframe(gcf);
    im_vectors = im_vectors.cdata;
    imwrite(im_vectors, fullfile(directory, 'images', ...
        ['piv_raw_', output_name, '.tif']), ...
        'writemode', 'append');
    
end
