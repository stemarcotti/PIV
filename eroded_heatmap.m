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
dilationSize = 4;       % [px]
erosionSize = 12;       % [px]
connectivityFill = 4;   % [px]
vector_spacing = 25;    % [px]
colour_min_val = 0;     % [um/min]

% read in the interpolated flow field
vfilt = load (fullfile ([d '/data'], ['piv_field_interpolated_', output_name, '.mat']));
vfilt = vfilt.vfilt;

%% ERODED HEATMAP %%

for jj = 1:length(vfilt)

    % read movie frames at t0 and t+1
    currentFrame = double(imread(fullfile(d, sprintf ...
        ('cb%d_m.tif', mt)),jj)) / 255;
    nextFrame = double(imread(fullfile(d, sprintf ...
        ('cb%d_m.tif', mt)),jj+1)) / 255;
    
    % read interpolated flow in x and y
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
    
    % read in movie without cell body (if available)
    file_name = [d, '/', sprintf('no_cb%d_m.tif', mt)];
    if exist(file_name, 'file') == 2
        no_cb_frame = double(imread(fullfile(d, sprintf ...
            ('no_cb%d_m.tif', mt)),jj)) / 255;
        
        % mask out cell body to vector field
        lim = logical(no_cb_frame);
        maskFlow.vx = maskFlow.vx .* lim;
        maskFlow.vy = maskFlow.vy .* lim;
        
    end
    
    % get the magnitude of the flow vectors
    vmag = sqrt(u.^2 + v.^2);
    vmag_no_cb = sqrt(maskFlow.vx.^2 + maskFlow.vy.^2);
    
    % plot
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

    % overlay black unit vectors onto flow field
    [x_str, y_str] = meshgrid(1:size(u,2), 1:size(u,1));

    quiver(x_str(1:vector_spacing:end,1:vector_spacing:end),...
        y_str(1:vector_spacing:end,1:vector_spacing:end),...
        maskFlow.vx(1:vector_spacing:end,1:vector_spacing:end)./vmag_no_cb(1:vector_spacing:end,1:vector_spacing:end),...
        maskFlow.vy(1:vector_spacing:end,1:vector_spacing:end)./vmag_no_cb(1:vector_spacing:end,1:vector_spacing:end),...
        0.3, 'k', ...
        'LineWidth', 1)
    
    % set figure window background colour to white
    set(gcf, 'InvertHardCopy', 'off');
    set(gcf, 'Color', [1 1 1]); 
    pause(1)
    hold off
    
    % save stack in folder [images]
    im_flow = getframe(gcf);
    im_flow = im_flow.cdata;

    imwrite(im_flow, fullfile(d, 'images', ...
        ['piv_interpolated_eroded_', output_name, '.tif']), ...
        'writemode', 'append');
    
    % this saves the movie at high resolution in individual tiff images
    % the save file path and the movie id needs to be manually added below
    
    %     % (change -dtiff to -dpng for less memory intensive images)
    %     print(gcf, '-dtiffn', '-r200', [d '/images/piv_interpolated_eroded_', output_name, '_frame' jj '_HR.tif']);
    %     pause(1)
    %     hold off
    
end
clear; close all