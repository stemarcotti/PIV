## PIV documentation

Particle Image Velocimetry (PIV) package developed in the Stramer Lab (King's College London, UK).

Tested on MATLAB v2018B. The Curve Fitting Toolbox is required.

This PIV code comes without a graphic user interface (GUI) and should be run in MATLAB as a script (open .m file and hit Run). It is used as follow.
See [here](https://www.sciencedirect.com/science/article/pii/S0092867415001816?via%3Dihub#app2) for references.

### 1. Image pre-processing
Segmented stacks of biological samples to be analysed should be pre-processed as follows:

1. Open stack in ImageJ
2. Split channels (e.g. green - actin, magenta - nucleus)
3. Save channel containing the entity to be measured by PIV (e.g. green - actin) as **[cb#\_m.tif]**, where cb stands for cell body, # is a progressive integer number, and m stands for moving
4. Save channel containing the entity used for tracking (e.g. magenta - nucleus) as **[n#_m.tif]**
5. If working with cells, segment out the cell body from the entity to be measured by PIV (e.g. green - actin) and save this one-channel stack as **[no_cb#\_m.tif]** (no cell body). This is needed for the **[eroded_heatmap.m]** script (see below).
6. Store **[cb#\_m.tif]**, **[n#_m.tif]** and **[no_cb#\_m.tif]** in a folder (e.g. _[cell1]_), where the PIV output is going to be saved

### 2. PIV
Clone this repository locally in your machine. Make sure the MATLAB folder in which you set your environment contains the following scripts: **[happy_piv.m]**, **[run_piv.m]**, **[create_retro_flow_image.m]**, **[detectObjectBw.m]**, **[generate_plot_normalized.m]** and, if needed, **[eroded_heatmap.m]**.

To run the PIV follow these instructions:

1. Run **[happy_piv.m]**
  - input requested to the user:
    + stack of the entity to be measured by PIV **[cb#\_m.tif]**
    + name for the output stamp to be appended to all saved output files (e.g. [output_name]: cell1)
    + PIV parameters:
      - _Source size [um]_: size of the region of interest (ROI) to be tracked in the next frame [um]
      - _Search size [um]_: size of the region to be used in the next frame for searching the ROI [um]
      - _Grid distance [um]_: grid spacing for ROI mapping [um]
      - _Correlation threshold [-]_
      - _Pixel size [um]_: image calibration of the loaded movie (length of a pixel in [um])
      - _Frame interval [s]_: time calibration of the loaded movie (frame interval in [s])
      - _Frame rate to be analysed [-]_: interval of frames of the loaded movie to be analysed (if 1: all the frames will be analysed, if 2: every other frame, and so on)
      - _Max number of frames to be analysed [frames]_ : maximum number of frames to be analysed (the default value displays the total available frames, do not change if wanting to analyse the whole movie)
    + Interpolation parameters (only prompted if interpolation is requested):
      - _Spatial kernel size [um]_: size of the spatial Gaussian kernel for interpolation
      - _Spatial kernel sigma [um]_: sigma (variance) of the spatial Gaussian kernel for interpolation
      - _Max flow velocity to be displayed in colourmap [um/min]_: upper limit of colour map in [um/min] for the heatmap
      - _Flow field arrow distance [um]_: distance of the arrow showing the vector field in the heatmap [um]
  - if parameters need to be optimised it is possible to show the frame correlation for a single frame, by hitting _Yes_ when prompted with the relevant question _"Do you want to test frame correlation to adjust the parameters?"_. This can be done multiple times, until the user doesn't reply _Yes_ to the following question _"Are you happy with the frame correlation test?"_.
  - Once the parameters are set, the PIV on all frames can be run by hitting _No_ when prompted with the question _"Do you want to test frame correlation to adjust the parameters?"_. The parameters are saved in the folder [parameters] which can be found where the loaded movie is stored (e.g. _[cell1]_).
  - The first phase of the PIV works out the raw vector fields in [um/min], which is visually displayed for each frame. The stack (**[piv_raw_(output_name).tif]**) and the numerical data (**[piv_field_raw_(output_name).mat]**) are saved in the folder [images] and [data] respectively. Both folders can be found where the loaded movie is stored (e.g. _[cell1]_).
  - The second phase of the PIV interpolates the raw vector fields to obtained the interpolated field heatmap in [um/min], which is visually displayed for each frame. This second part is optional, and will be run only if the user answer _Yes_ when prompted with the question _"Do you want to interpolate the vector field?"_. The spatial kernel parameters are chosen by the user, while the temporal kernel interpolates within 5 frames as default. The stack (**[piv_interpolated_(output_name).tif]**) and the numerical data (**[piv_field_interpolated_(output_name).mat]**) are saved in the folder [images] and [data] respectively.  

2. More refined interpolated heatmaps can be obtained by running the **[eroded_heatmap.m]** script
  - input requested to the user:
    + folder containing **[cb#\_m.tif]**, **[no_cb#\_m.tif]** and PIV output (e.g. _[cell1]_). If working with cell images for which the **[no_cb#\_m.tif]** is available, the final heatmaps will not show overlayed vectors for this area. If the file **[no_cb#\_m.tif]** does not exist, the heatmap will display the original entity in full.
    + name for the output stamp to be appended to all saved output files (e.g. [output_name]: cell1); need to be the same assigned when running **[happy_piv.m]**!
    + the movie ID (# in **[cb#\_m.tif]** and **[no_cb#\_m.tif]**)
  - this script returns the refined stack (**[piv_field_interpolated_eroded_(output_name).tif]**), containing the heatmaps of the interpolated PIV field overlayed with the flow vectors (excluded the cell body if applicable).
