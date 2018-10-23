%% INPUT %%

warning off

% load cell movie (user selects folder and .tif file)
uiwait(msgbox('Load cell movie'));
[file,directory] = uigetfile('.tif');

% ask the user for an ouput stamp
prompt = {'Provide a name for the output files'};
title = 'Output name';
dims = [1 35];
user_answer = inputdlg(prompt,title,dims);
output_name = (user_answer{1,1});

% number of frames of the selected movie
nt = length(imfinfo(fullfile(directory,file)));

% create output folders if not already there
if ~exist(fullfile(directory, 'data'))
    mkdir(fullfile(directory, 'data'))
end

if ~exist(fullfile(directory, 'images'))
    mkdir(fullfile(directory, 'images'));
end

if ~exist(fullfile(directory, 'parameters'))
    mkdir(fullfile(directory, 'parameters'));
end


%% RUN PIV and obtain raw output %%

% initialise while loop
happy_test = 0;

% this while loop will run until the user is happy with the correlation
% frame test (until [happy_test] = 1)
while happy_test == 0

    % user input parameters for PIV
    prompt = {'Source size [um]', ...
        'Search size [um]',...
        'Grid distance [um]',...
        'Correlation threshold [-]',...
        'Pixel size [um]',...
        'Frame interval [s]',...
        'Frame rate to be analysed [-]',...
        'Max number of frames to be analysed (default displays total available frames) [frames]'};
    title = 'PIV parameters';
    dims = [1 35];
    definput = {'1.2','2.0','0.8','0.5', '0.1', '5', '1', num2str(nt)};
    user_answer = inputdlg(prompt,title,dims,definput);

    source_size_user = str2double(user_answer{1,1});              % source size [um]
    search_size_user = str2double(user_answer{2,1});              % search size [um]
    grid_distance_user = str2double(user_answer{3,1});            % grid distance [um]
    params.correlation_threshold = str2double(user_answer{4,1});	% correlation threshold [-]
    params.mu2px = str2double(user_answer{5,1});                    % pixel size [um]
    params.max_frame = str2double(user_answer{8,1});                % max number of frames to be analysed

    params.frame_rate = str2double(user_answer{7,1});               % frame rate to be analysed

    params.recording_speed = str2double(user_answer{6,1}) * params.frame_rate;  % recording speed (frame interval [s])
    fps = 1 / params.recording_speed;	% frames per second [frame/s]

    params.source_size = source_size_user/params.mu2px;  % [px]
    params.search_size = search_size_user/params.mu2px;  % [px]
    params.grid_distance = grid_distance_user/params.mu2px;  %[px]

    % ask the user if they want to run the test frame correlation on the
    % first frame
    run_test_q = questdlg('Do you want to test frame correlation to adjust the parameters?', ...
        'Run test',...
        'Yes', 'No', 'No');

    if strcmp(run_test_q, 'Yes')
        run_test = 1;	% if yes: [run_test] is assigned to the value 1
        k = 1;          % and k = 1 means the PIV only runs on frame 1
        run_piv(directory, file, k, params, happy_test, output_name);
    else
        run_test = 0;   % if no: [run_test] equals 0
        max_k = params.max_frame-1;   % and the PIV is run on all the frames of the movie

        happy_test = 1; % this will allow to exit the while loop

        vraw_index = 1; % this is needed to save the [vraw] output

        frame_counter1 = 1;
        for k = 1:params.frame_rate:max_k

            [vraw(vraw_index).x, vraw(vraw_index).y, vraw(vraw_index).vx, vraw(vraw_index).vy, vraw(vraw_index).cc] = run_piv(directory, file, k, params, happy_test, output_name);
            vraw_index = vraw_index+1;

            fprintf('Running PIV: frame %d/%d \n', frame_counter1, ceil(max_k/params.frame_rate))
            frame_counter1 = frame_counter1 + 1;

        end
    end

    uiwait(msgbox('Press OK to continue')); % wait user OK to close figure
    close all

    % if the user is running the test: ask for feedback
    if run_test == 1
        happy_test_q = questdlg('Are you happy with the frame correlation test?', ...
            'Happy?',...
            'Yes', 'No', 'Yes');

        if strcmp(happy_test_q, 'Yes')
            happy_test = 1;	% if yes: [happy_test] is assigned to the value 1 (exit while loop)
            max_k = params.max_frame-1;   % and the PIV for all the frames is computed

            vraw_index = 1;
            frame_counter1 = 1;
            for k = 1:params.frame_rate:max_k

                [vraw(vraw_index).x, vraw(vraw_index).y, vraw(vraw_index).vx, vraw(vraw_index).vy, vraw(vraw_index).cc] = run_piv(directory, file, k, params, happy_test, output_name);
                vraw_index = vraw_index+1;

                fprintf('Running PIV: frame %d/%d \n', frame_counter1, ceil(max_k/params.frame_rate))
                frame_counter1 = frame_counter1 + 1;
            end

            uiwait(msgbox('Press OK to continue'));	% wait user OK to close figure
            close all

        else
            happy_test = 0; % if no: run test again (don't exit the while loop)
        end
    end
end

%% raw OUTPUT %%

% save raw PIV field in [data]
save(fullfile([directory '/data'], ...
    ['piv_field_raw_', output_name, '.mat']), ...
    'vraw');

%% FILTER and INTERPOLATE %%
fprintf('\n')

% ask the user if they want to run the interpolation
interp_q = questdlg('Do you want to interpolate the vector field?', ...
    'Interpolation',...
    'Yes', 'No', 'Yes');

if strcmp(interp_q, 'Yes')

    % user input parameters for interpolation
    prompt = {'Spatial kernel size [um]', ...
        'Spatial kernel sigma [um]',...
        'Max flow velocity to be displayed in colourmap [um/min]',...
        'Flow field arrow distance [um]'};
    title = 'Interpolation parameters';
    dims = [1 35];
    definput = {'5', '1', '12', '2'};
    user_answer = inputdlg(prompt,title,dims,definput);
    params.k_size_user = str2double(user_answer{1,1});  % spatial kernel size [um]
    params.k_sigma_user = str2double(user_answer{2,1}); % spatial kernel sigma [um]
    params.max_flow_vel = str2double(user_answer{3,1});	% max flow velocity [um/min]
    params.flow_field_arrow_distance = str2double(user_answer{4,1})/params.mu2px;	% arrow spatial resolution [px]

    % change unit to parameters
    k_size = 2 * (ceil(0.5 * params.k_size_user / params.mu2px));	% spatial kernel size [px]
    k_sigma = 2 * (ceil(0.5 * params.k_sigma_user / params.mu2px));	% spatial kernel sigma [px]
    k_size_temp = 4;   % temporal kernel size [frames] (-1 is added because it sums up 1 afterwards)
    k_sigma_temp = 2;       % temporal kernel sigma [frames]

    % switch from pixel per image to mu/min
    for i = 1:length(vraw)
        list(:, 1) = vraw(i).x; % [px]
        list(:, 2) = vraw(i).y; % [px]
        list(:, 3) = vraw(i).vx * params.mu2px * fps * 60;  % [um/min]
        list(:, 4) = vraw(i).vy * params.mu2px * fps * 60;  % [um/min]
        list(:, 5) = vraw(i).cc; % [-]
        list_cell{i} = list;
        clear list
    end

    % get the max and min x and y values of the flow vectors startpoints for
    % the full image series (means you get the values closest to the image
    % boarder in all the time series)
    t = max(size(list_cell));
    for j = 1:t
        list = list_cell{j};
        y(j) = max(list(:, 2));
        x(j) = max(list(:, 1));
        y_min(j) = min(list(:, 2));
        x_min(j) = min(list(:, 1));
    end

    % to prevent problems with flow values too close to the upper or left image
    % boarder (so that the filter kernel would reach outside the image and
    % negative array indices are not defined, in the positive direction there
    % is no problem, the code just fills up the missing values with zeros),
    % check if the min x,y values of the flow points in list_cell are within
    % kernel/2. If not, shift the startpoints of the flow arrows, and all the
    % images in the series. All the shift info is stored in the x_shift, and
    % y_shift variables.
    if (min(x_min) <= k_size/2)
        x_shift = ceil(k_size/2 - min(x_min) + 1);
    else
        x_shift=0;
    end

    if (min(y_min) <= k_size/2)
        y_shift = ceil(k_size/2 - min(y_min) + 1);
    else
        y_shift=0;
    end

    % shift the retroflow values for the full image series.
    for j = 1:t
        list = list_cell{j};
        list(:,2) = list(:,2) + y_shift;
        list(:,1) = list(:,1) + x_shift;
        list_cell{j} = list;
    end

    % define the k_size_temp+1 (e.g 5) layer thick convolution array block
    % that is moved through the full image series in the analysis.
    x = max(x);
    y = max(y);
    im_x(1:y + k_size + 1 + y_shift, 2:x + k_size + 1 + x_shift, ...
        1:k_size_temp + 1) = 0;
    im_divider = im_x;
    im_y = im_x;

    % generate the convolution kernel, gaussian kernel, circular,
    % sigma of k_size/6.
    [k_x, k_y] = meshgrid(-k_size/2:k_size/2, -k_size/2:k_size/2);
    kernel_2D = exp(-(k_x.^2 + k_y.^2) / (2 * (k_sigma)^2));    % [um]
    for i = 1:k_size_temp/2 + 1
        kernel_3D(:,:,i) = kernel_2D * exp(-(k_size_temp/2 + 1 - i)^2 / ...
            (k_sigma_temp)^2);   % [um]
        kernel_3D(:,:,k_size_temp - i + 2) = kernel_3D(:,:,i);
    end
    min_val = max(min(kernel_3D(:,:,k_size_temp/2 + 1)));
    kernel_3D(kernel_3D < min_val) = 0;

    startframe = 1;
    endframe = max(size(list_cell));

    % 3D convolution
    figure()
    axes = gca;
    frame_counter2 = 1;
    for j = startframe:endframe
        list = list_cell{j};

        fprintf('Interpolating: frame %d/%d \n', frame_counter2, endframe)
        frame_counter2 = frame_counter2 + 1;

        for i = 1:max(size(list))
            x_c = list(i,1);
            y_c = list(i,2);
            dx = list(i,3);
            dy = list(i,4);

            if(list(i,5) > params.correlation_threshold)
                x_c = list(i,1);
                y_c = list(i,2);

                % build the convolution matrix and the normalization matrices
                im_x(y_c - (k_size/2):y_c + (k_size/2), ...
                    x_c - (k_size/2):x_c + (k_size/2), :) = ...
                    list(i,5) * list(i,3) * kernel_3D + ...
                    im_x(y_c - (k_size/2):y_c + (k_size/2), ...
                    x_c - (k_size/2):x_c + (k_size/2), :); % [um2/s]

                im_y(y_c - (k_size/2):y_c + (k_size/2), ...
                    x_c - (k_size/2):x_c + (k_size/2), :) = ...
                    list(i,5) * list(i,4) * kernel_3D + ...
                    im_y(y_c - (k_size/2):y_c + (k_size/2), ...
                    x_c - (k_size/2):x_c + (k_size/2), :); % [um2/s]

                im_divider(y_c - (k_size/2):y_c + (k_size/2), ...
                    x_c - (k_size/2):x_c + (k_size/2), :) = ...
                    list(i,5) * kernel_3D + ...
                    im_divider(y_c - (k_size/2):y_c + (k_size/2), ...
                    x_c - (k_size/2):x_c + (k_size/2), :); % [um]
            end
        end

        im_x_act = single(im_x(:,:,1) ./ im_divider(:,:,1)); % [um/s]
        im_y_act = single(im_y(:,:,1) ./ im_divider(:,:,1));

        % create image and save final values (starts at frame k_size_temp/2 + 1
        % because of time filtering)
        if j > startframe-1+k_size_temp/2

            i = j - k_size_temp / 2;

            % create the image, plus get the filtered final results
            % load current frame
            im = imread(fullfile(directory, file), i);
            im = double(im) / 255;

            % detect cell outline for current frame
            edge_alg = 'Canny';
            im_edge = edge(im, 'Canny');
            im_edge = imdilate(im_edge, ...
                strel('disk', 4));
            im_edge = imfill(im_edge, 4, 'holes');
            im_edge = imerode(im_edge, ...
                strel('disk', 4));

            % create flow heatmaps
            [vfilt(i).vx, vfilt(i).vy] = create_retro_flow_image(axes, ...
                directory, im, im_edge, im_x_act, im_y_act, x_shift, ...
                y_shift, params.mu2px, params.max_flow_vel, ...
                params.flow_field_arrow_distance, output_name);
        end

        % delete first time layer of convolution array block and define next one
        im_x(:,:,1) = [];
        im_y(:,:,1) = [];
        im_divider(:,:,1) = [];
        im_x(:,:,k_size_temp + 1) = 0;
        im_y(:,:,k_size_temp + 1) = 0;
        im_divider(:,:,k_size_temp + 1) = 0;

        % create the last image save final values (has to be done differently because of time filtering)
        if j == max(size(list_cell)) % this means the last 2 frames run again but its easy this way

            for n = 1:k_size_temp/2
                i = j - k_size_temp / 2 + n;
                im_x_act = single(im_x(:,:,n) ./ im_divider(:,:,n));
                im_y_act = single(im_y(:,:,n) ./ im_divider(:,:,n));

                %  create the image, plus get the filtered final results
                [vfilt(i).vx, vfilt(i).vy] = create_retro_flow_image(axes, ...
                    directory, im, im_edge, im_x_act, im_y_act, x_shift, ...
                    y_shift, params.mu2px, params.max_flow_vel, ...
                    params.flow_field_arrow_distance, output_name);
            end
        end
    end

    %% interpolated OUTPUT %%

    % save interpolated flow data to .mat file
    save(fullfile(directory, 'data', ...
        ['piv_field_interpolated_', output_name,'.mat']), ...
        'vfilt');


    % save PIV parameters in [parameters]
    save(fullfile([directory '/parameters'], ...
        ['piv_parameters_', output_name, '.mat']), ...
        'params');

end

clear; close all;
