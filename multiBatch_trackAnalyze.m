%% Batch Analyze Tracks from a Multiple Conditions

% Instructions
% Create a folder with the name you want to use as the condition name
% Put all the linkedTracks.mat files for that condition in the folder
% Put all the condition folders in a single folder, this is what you will
% direct the program to

% Notes
% This code does not interpolate/tracks are not linked

% Dependencies
% fix_tracks_allfields
% get_mat_from_tracks
% getSpeedAnalysis

%% Find the files for each condition

% Clear pre-existing data
clear

% Save the current path
analysisPath = pwd;

% Let user select directory
dataPath = uigetdir('Select directory');
cd(dataPath)
folders = dir();
folders = folders([folders.isdir]);  % Get only directories
folderNames = {folders(~ismember({folders.name}, {'.', '..'})).name};

%% Run the analysis for each folder

for i = 1:size(folderNames, 2)
    currentFolder = folderNames{i};
    miniTrackAnalyze(dataPath, currentFolder, analysisPath);
end 

