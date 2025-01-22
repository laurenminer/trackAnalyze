%% Batch Analyze Tracks from a Single Condition

% Instructions
% Create a folder with the name you want to use as the condition name
% Put all the linkedTracks.mat files for that condition in the folder

% Notes
% This code does not interpolate/tracks are not linked

% Dependencies
% fix_tracks_allfields
% get_mat_from_tracks
% getSpeedAnalysis

%% Find all the linkedTracks Files

% Clear pre-existing data
clear

% Save the current path
analysisPath = pwd;

% Let user select directory
dirPath = uigetdir('Select directory');
pathParts = split(dirPath, filesep);
conditionName = pathParts{end};

% Get all .mat files ending with linkedTracks.mat
filePattern = fullfile(dirPath, '*linkedTracks.mat');
files = dir(filePattern);
fileNames = string({files.name});

%% Extract, Concatenate, and Exportt= Data

% Go to the place where the data is 
cd(dirPath);

% Housekeeping
nFiles = length(fileNames);
unlinkedData = struct();

% fieldsToSave is a list of fields present in the .unlinkedData.mat data files to directly copy to chxData 
fieldsToSave = ["Frames","Time","Speed","AngSpeed","SmoothX","SmoothY","Curvature"]; % 
for f = 1:nFiles % each row of unlinkedData corresponds to a video file
    name = fileNames(f);

    % Load data. If it is a linkedTracks file you need to create and save
    % the unlinkedData file first.
    if contains(name,"linkedTracks.mat")
        load(name);
        unlinkedData = linkedTracks; % create unlinkedData.mat
        % THIS IS THE LINE TO EDIT FOR CHANGING HOW TRACKS CONNECT:
        unlinkedData = fix_tracks_allfields(unlinkedData); % remove random unwanted data 
        storeName = strrep(name,".linkedTracks.mat","");
    elseif contains(name,"linkedTracks.mat")
        load(name);
        storeName = strrep(name,".Tracks.mat","");
    else
        error("Must be a .linkedTracks.mat or .linkTracks.mat")
    end

    unlinkedData(f).Name = storeName;
    
    for field = fieldsToSave % loop through fields and store the data in time-aligned nTracks x nFrames matrices
        unlinkedData(f).(field) = single(get_mat_from_tracks(linkedTracks,field));
    end
end 

% Save the concantenated data file (with a new name if it already exists)
savename = strcat(dirPath, '\', conditionName,'.unlinkedData.mat');
counter = 1;
while true
    if exist(savename,"file") % don't overwrite
        savename = strcat(dirPath, '\', conditionName,'_version',num2str(counter),'.unlinkedData.mat');
        counter = counter+1;
    else
        try
            save(savename,'unlinkedData')
        catch
            warning('Error saving as a v7 file, probably due to large size. Trying v7.3.')
            save(savename,'unlinkedData','-v7.3')
        end
        break;
    end
end

% Go back to the analysis folder
cd(analysisPath)

%% Housekeeping

pathParts = split(savename, filesep);
unlinkedDataName = pathParts{end};

clearvars -except analysisPath dirPath conditionName unlinkedDataName

%% Extract Speed, Roaming, and Dwelling Data

getSpeedAnalysis(unlinkedDataName)