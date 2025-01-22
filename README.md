## The trackAnalyze Package processes a set of linkedTracks.mat files and reports back analyzed speed and roaming/dwelling state information.

# Prepare your data
To prepare your data for analysis, put all the linkedTracks files for a given condtion in a folder with nothing else. The name of this folder will be the condition name displayed on graphs. If you have single condition, you can use the singleBatch_trackAnalyze script. Navigate to that folder when you are prompted.

If you have multiple conditions you can put all these condition folders in a single folder (put nothing else in this folder) and use the multiBatch_trackAnalyze.script. Navigate to the folder that has all the conditions (not one of the subfolders) when prompted. multiBatch uses the same settings to analyze all conditions.

# What you get in return

In the folder for your condition you will get back 3 matlab variables:
1. In unlinkedData, each row is a video/linkedTracks file and the columns have the raw time, speed, angular speed, and curvature data from the original linkedTracks file in each column.
2. In binnedData, each row is a video/linkedTracks filer but the relevant data has been binned into a set interval.
3. In processedData, all the binned data has been grouped together. You can see this data as well as the settings used for key variables and the fraction of time spent roaming/dwelling calculated by different methods (for plotting elsewhere or recordkeeping).

You will get back a number of plots all saved both as matlab figures and .jpeg files. See below for more details.

All data and figures are saved in the folder for a condition.

# Customizations
All the parameters you would want to change depending on your experiment are in getSpeedAnalysis under the Plot Settings and Variables/Housekeeping sections.

1. Plotting
   There are many options for what to plot. In getSpeedAnalysis at the top of the code you can choose to proceed plotting and saving a graph by setting the corresponding variable equal to 1, or not by setting it to 0. These are the plots and the corresponding variable to change in getSpeedAnalysis
   - Scatter Plot of Speed vs Time: speed_time_scatter
   - Heatmap of Speed vs Time: speed_time_heatmap
   - Scatter Plot of Angular Speed vs Time: ang_time_scatter
   - Heatmap of Angular Speed vs Time: ang_time_heatmap
   - Scatter Plot of Speed vs Angular Speed: speed_ang_scatter
   - Heatmap of Speed vs Angular Speed: speed_ang_heatmap
   - Scatter Plot of Speed vs Curvature: speed_curve_scatter
   - Heatmap of Speed vs Curvature: speed_curve_heatmap
   - XY plot of Speed vs Time where each line is a plate (single linkedTracksFile): speed_time_plate
   - Stacked Bar of Fraction Spent Roaming vs Dwelling, Calulated by a HMM, Roaming and Dwelling Observations informed by the relationship between Speed and Angular Speed (Flavell 2013): roam_dwell_speed_ang
   - Stacked Bar of Fraction Spent Roaming vs Dwelling, Calulated by a HMM, Roaming and Dwelling Observations informed by a speed threshold: roam_dwell_only_speed
   - Stacked Bar of Fraction Spent Roaming vs Dwelling, Calulated by a HMM, Roaming and Dwelling Observations informed by the relationship between Speed and Curvature (Ben Arous 2009): roam_dwell_speed_curve
  
    All the figures will pop up as they are produced. If you want them to stay up, set closeThem = 0. If you want them to close, set closeThem = 1. They are saved either way.

2. Classifying Roaming vs Dwelling Observations (setPoint, speedPoint, and curvePoint variables)
   The animal's behavior in any given time bin can be classified in three ways (at least for this code):
   1. In Flavell 2013, speed was plotted against angular velocity and a line (y = x/setPoint) was drawn. Observations above the line are classified as roaming. Observations below the line are classified as dwelling. To adjust the line, change the setPoint variable (getSpeedAnalysis line 40). The line is automatically plotted in the Speed vs Angular Speed Heatmap. In the paper, the inital observations used a setPoint of 450, but this could change from experiment to experiment and a range of values from 270 to 550 were used.
   2. You could set a threshold where all observations above a certain speed are classified as roaming, all below are dwelling. To change this threshold, change the speedPoint Variable ((getSpeedAnalysis line 41). The line is automatically plotted in the Speed vs Time Heatmap.
   3. In Ben Arous 2009, speed was speed was plotted against curvature and a line (y = x/curvePoint) was drawn. Observations above the line are classified as roaming. Observations below the line are classified as dwelling. To adjust the line, change the setPoint variable (getSpeedAnalysis line 42). The line is automatically plotted in the Speed vs Curvature Heatmap. In the paper, the curvePoint used looks to be approximately 200. This was used for all analysis.

3. Other Parameters You Need to Look At
   1. frameRate is the frame rate of your camera. Our camera is 3 fps, this is the preset value.
   2. binSec is how many seconds you want to group in a bin (the data in this bin gets averaged together to make one point). Both Flavell 2013 and Ben Arous 2009 use 10 second bins, this is the preset value. Bins are pre-established for all videos based on linear time (not based on when you start having non-nan data in a track, see iv below if you don't know what that means).
   3. totalFrames is how long you expect your video to be. If an input video/linkedTracks file is longer or shorter than that, the code will cut off extra data or pad the data with NaN values (not plotted) so you don't get errors. My videos are 1 hr/10800 frames long so this is the preset value.
   4. nanLimit. Sometimes the tracker loses track of a worm or worms overlap. This data is not interpolated and instead shows up as NaN values. The nanLimit is how many nan frames you want to allow in a bin and still use that bin. If the code does not use that bin, it will replace whatever it would have calulated in that bin with NaN for all calculations. If you have the same frame rate (3) and bin size (10s) as me, then there are 30 frames in a bin. I only allow 3 of those values to be NaN (90% of the data is there), so 3 is the preset. This was an arbitrary decision.
   
