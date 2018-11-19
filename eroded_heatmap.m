%% INPUT %%
% get the file directory
uiwait(msgbox('Load cell movie folder'));
d = uigetdir('');
listing = dir (fullfile (d, 'cb*.tif'));
numFiles = length (listing);

% ask the user for an ouput stamp
prompt = {'Provide a name for the output files',...
    'Movie ID (n) if file format is cb_(n)_m.tif',...
    'Max flow velocity to be displayed in colourmap [um/min]'};
title = 'Parameters';
dims = [1 35];
user_answer = inputdlg(prompt,title,dims);
output_name = (user_answer{1,1});
mt = str2double(user_answer{2,1});
colour_max_val = str2double(user_answer{3,1}); % [um/min]

% parameters
dilationSize = 4;  % [px]
erosionSize = 12;  % [px]
connectivityFill = 4;  % [px]
div_filt_size = 3 * dilationSize;
is_erode = 1;
outline_dilation_size = 4;  % [px]
outline_erosion_size = 1;   % [px]
vector_spacing = 25;  % [px]
colour_min_val = 0;   % [um/min]

% read in the interpolated flow field
vfilt = load (fullfile ([d '/data'], ['piv_field_interpolated_', output_name, '.mat']));
vfilt = vfilt.vfilt;

%% ERODE HEATMAP %%
for jj = 1:length(vfilt)

    % read movie frames at t0 and t+1
    currentFrame = double(imread(fullfile(d, sprintf ...
        ('cb%d_m.tif', mt)),jj)) / 255;
    nextFrame = double(imread(fullfile(d, sprintf ...
        ('cb%d_m.tif', mt)),jj+1)) / 255;

    u = vfilt(jj).vx;
    v = vfilt(jj).vy;

    % produce masks of movie at t0 and t+1
    cellOutline1 = detectObjectBw(currentFrame, dilationSize, erosionSize, connectivityFill);
    cellOutline2 = detectObjectBw(nextFrame, dilationSize, erosionSize, connectivityFill);

    % get the intersection of two consecutuve masks
    cellOutline = cellOutline1 .* cellOutline2;

    % mask off the flow field (it will remove 'blobby' edges)
    maskFlow.vx = u .* cellOutline;
    maskFlow.vy = v .* cellOutline;
    % read in movie without cell body
    file_name = [d, '/', sprintf('no_cb%d_m.tif', mt)];
    if exist(file_name, 'file') == 2
        no_cb_frame = double(imread(fullfile(d, sprintf ...
            ('no_cb%d_m.tif', mt)),jj)) / 255;

        % convert the frame without the cell body to a mask
        lim = logical(no_cb_frame);
            maskFlow.vx = maskFlow.vx .* lim;
    maskFlow.vy = maskFlow.vy .* lim;

    end




    % get the magnitude of the flow vectors
    vmag = sqrt(u.^2 + v.^2);

    % get the magnitude of the flow vectors
    vmag_n = sqrt(maskFlow.vx.^2 + maskFlow.vy.^2);

    % find cell outline
    im = currentFrame .* cellOutline;
    imbw = edge(im, 'canny');
    imbw = imdilate(imbw, strel('disk', outline_dilation_size));
    imbw = imfill(imbw, 'holes');
    imbw = imerode(imbw, strel('disk', outline_dilation_size));
    imbw = double(imbw);

    % erode cell outline to get rid of anomalous values around the cell edge
    if is_erode == 1
        im_filt = imerode(imbw, strel('disk', outline_erosion_size));
        im_filt = double(im_filt);
    end

    vectors = maskFlow;
    if ~isempty(vectors) && is_erode == 1
        vx = vectors.vx .* im_filt;
        vy = vectors.vy .* im_filt;
        clear vectors
    else
        vx = vectors.vx;
        vy = vectors.vy;
        clear vectors
    end

    % display image
    data = vmag .* cellOutline;

    h = imshow(data, []);

    colormap('jet');
    caxis([colour_min_val, colour_max_val])
    c = colorbar;
    c.Label.FontSize = 14;
    c.Label.String = 'Velocity [um/min]';
    hold on

    % black background
    data(data == 0) = NaN;
    set(h, 'AlphaData', ~isnan(data))
    axis on;
    set(gca, 'XColor', 'none', 'yColor', 'none', 'xtick', [], 'ytick', [], 'Color', 'k')
    hold on

    [x_str, y_str] = meshgrid(1:size(u,2), 1:size(u,1));

    % overlay black unit vectors onto flow field
    quiver(x_str(1:vector_spacing:end,1:vector_spacing:end),...
        y_str(1:vector_spacing:end,1:vector_spacing:end),...
        maskFlow.vx(1:vector_spacing:end,1:vector_spacing:end)./vmag_n(1:vector_spacing:end,1:vector_spacing:end),...
        maskFlow.vy(1:vector_spacing:end,1:vector_spacing:end)./vmag_n(1:vector_spacing:end,1:vector_spacing:end),...
        0.3,'w', 'LineWidth',1)

    set(gcf, 'InvertHardCopy', 'off');
    set(gcf, 'Color', [1 1 1]); %setting figure window background color white
    pause(1)
    hold off

    % this saves the movie at high resolution in individual tiff images
    % the save file path and the movie id needs to be manually added below

%     % (change -dtiff to -dpng for less memory intensive images)
%     print(gcf, '-dtiffn', '-r200', [d sprintf('/Movie_%d_Flow_Field%d',mt, jj)]);
%     pause(1)
%     hold off

    % save stack in folder [images]
    im_flow = getframe(gcf);
    im_flow = im_flow.cdata;

    imwrite(im_flow, fullfile(d, 'images', ...
        ['piv_interpolated_eroded_', output_name, '.tif']), ...
        'writemode', 'append');
%     close all;

end
clear; close all
