function [x_flow, y_flow] = create_retro_flow_image(fig_handle, rf_path, im, ...
    bin_im, x_flow, y_flow, x_shift, y_shift, mue2pix_ratio, max_r_flow, ...
    flow_field_arrow_distance, output_name)
% Create heatmap with flow vectors.
%
% Input:
%   fig_handle                = handle for the figure where to plot the
%                               heatmaps
%   rf_path                   = path to the root folder where to save
%                               heatmaps. The code will automatically
%                               save a .tif file with the heatmaps in
%                               rf_path/images/
%   im                        = image with the current frame.
%   bin_im                    = image with the cell outline in the current
%                               frame.
%   x_flow                    = [M x N] matrix with the x-component of the
%                               flow field as produced by interpolation.
%   y_flow                    = [M x N] matrix with the y-component of the
%                               flow field as produced by interpolation.
%   x_shift                   = scalar with how much to shift the
%                               flow field in x direction. If running
%                               create_retro_flow_image after interpolation,
%                               shift should be 0.
%   y_shift                   = scalar with how much to shift the
%                               flow filed in y direction. If running
%                               create_retro_flow_image after interpolation,
%                               shift should be 0.
%   mue2pix_ratio             = microns per pixel ratio
%   max_r_flow                = maximum value of the flow. Used to scale
%                               displayed values between [0, max_r_flow].
%   flow_field_arrow_distance = distance in pixels at which to plot flow
%                               field vectors.
%   token                     = char specifying the type of flow data.
%                               Allowed values: 'm' = moving cell;
%                               's' = stationary cell.
%   running_tag               = string to tag final .tif file. Use the
%                               following format: '_HH_MM_DD_MM_YYYY'
%   output_name               = output name provided by the user
%
% Output:
%   x_flow                    = [M x N] matrix with the scaled flow field.
%   y_flow                    = [M x N] matrix with the scaled flow field.
%
% Note: to produce heatmaps with the flow field for an image sequence,
% create_retro_flow_image should be run in a loop for every frame.

warning off MATLAB:divideByZero
axes(fig_handle);

% Get the image size.
[y_max, x_max] = size(im);

% Take care of the x_shift, y_shift;
im_int(1:y_shift + y_max, x_shift + x_max) = 0;
im_int(y_shift + 1:y_shift + y_max, x_shift + 1:x_shift + x_max) = im;
im = im_int;

% Take care of the image shift
im_int(1:y_shift + y_max, x_shift + x_max) = 0;
im_int(y_shift + 1:y_shift + y_max, x_shift + 1:x_shift + x_max) = bin_im;
r_im_erode = im_int;

[y_max, x_max] = size(im);
max_size = max(size(im));

x_flow = double(x_flow);
y_flow = double(y_flow);
[y_max_field, x_max_field] = size(x_flow);

% Make sure image and the flow data have the same size; if not reshape
% the flow data accordingly.
if size(x_flow, 1) ~= size(im, 1) || size(x_flow, 2) ~= size(im, 2)
  intermediate(1:y_max, 1:x_max) = NaN;
  intermediate(1:min([y_max, y_max_field]), ...
      1:min([x_max, x_max_field])) = ...
      x_flow(1:min([y_max, y_max_field]), 1:min([x_max, x_max_field]));
  x_flow = intermediate;

  intermediate(1:y_max, 1:x_max) = NaN;
  intermediate(1:min([y_max, y_max_field]), ...
      1:min([x_max, x_max_field])) = ...
      y_flow(1:min([y_max, y_max_field]), ...
      1:min([x_max, x_max_field]));
  y_flow = intermediate;
end
clear('intermediate');

% Use the shape, to cut of whatever of the retro flow maps is outside
% the real shape.
x_flow = r_im_erode .* x_flow;
y_flow = r_im_erode .* y_flow;
x_flow(isnan(x_flow)) = 0;
y_flow(isnan(y_flow)) = 0;
x_flow_i = x_flow(x_shift + 1:end, y_shift + 1:end);
y_flow_i = y_flow(x_shift + 1:end, y_shift + 1:end);
x_flow = x_flow_i;
y_flow = y_flow_i;
[~, r_flow] = cart2pol(x_flow, y_flow);

% Drawing flow heatmaps.
% Make the color coded flow.
r_flow(isnan(r_flow)) = 0;
imshow(r_flow);
colormap('jet')
set(gca, 'CLim', [0, max_r_flow])
hold on

% Draw the flow field arrows.
plot_data = generate_plot_normalized(x_flow, y_flow, ...
    flow_field_arrow_distance);
quiver(plot_data(:,1), plot_data(:,2), plot_data(:,3), plot_data(:,4), ...
    0.3, 'k');

% Draw a 10 um scale bar
fill([x_max - x_shift - 5 - 10 / mue2pix_ratio, ...
    x_max - x_shift - 5 - 10 / mue2pix_ratio, ...
    x_max - 5 - x_shift, x_max - 5], [5, 10, 10, 5], 'w');

% Save the images, the display_edge_array, and of course the
% final_results_out.

h = figure('visible','off');
imshow(r_flow);
colormap('jet')
set(gca, 'CLim', [0, max_r_flow])
hold on

% Draw the flow field arrows.
plot_data = generate_plot_normalized(x_flow, y_flow, ...
    flow_field_arrow_distance);
quiver(plot_data(:,1), plot_data(:,2), plot_data(:,3), plot_data(:,4), ...
    0.3, 'k');

% Draw a 10 um scale bar
fill([x_max - x_shift - 5 - 10 / mue2pix_ratio, ...
    x_max - x_shift - 5 - 10 / mue2pix_ratio, ...
    x_max - 5 - x_shift, x_max - 5], [5, 10, 10, 5], 'w');

im_heat = getframe(h);
im_heat = im_heat.cdata;

imwrite(im_heat, fullfile(rf_path, 'images', ...
    ['piv_interpolated_', output_name, '.tif']), ...
    'writemode', 'append');

close(h);
