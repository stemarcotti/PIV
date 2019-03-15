function BW2 = detectObjectBw(im, dilationSize, erosionSize, connectivityFill)
    
    % return largest object only
    
    im = double(im) / 255;
    bw = edge(im, 'Canny');
    
    BW2d = imdilate(bw, strel('disk', dilationSize));
    BW2f = imfill(BW2d, connectivityFill, 'holes');
    BW2 = imerode(BW2f, strel('disk', erosionSize));
    
    [L,~] = bwlabel(BW2);
    
    stats = regionprops(L, 'Area');
    allArea = [stats.Area];
    area_largest_obj = max(allArea(:));
    
    BW2 = bwareaopen(BW2,area_largest_obj);

end
