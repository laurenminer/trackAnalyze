function getSpeedAnalysis(dirPath, unlinkedDataName, conditionName)


% Calculates the speed in 10 second intervals
% Intervals are pre-set as consecutive blocks of 30 frames
% Calculation is not performed if there are more than 3 frames missing

%% Load Data
cd(dirPath);
load(unlinkedDataName)

%% Plot settings
% 0 is do not plot
% 1 is plot
speed_time_scatter = 0; 
speed_time_heatmap = 0;
ang_time_scatter = 0;
ang_time_heatmap = 0;
speed_ang_scatter = 0;
speed_ang_heatmap = 1;
speed_curve_scatter = 0;
speed_curve_heatmap = 0;
speed_time_plate = 0;
roam_dwell_speed_ang = 0;
roam_dwell_only_speed = 0;
roam_dwell_speed_curve= 0;

% Do you want the plots to pop up and stay up or close
% 0 = keep them up
% 1 = close them
closeThem = 1;

%% Variables/Housekeeping

% Hard coded variables determined by experimental set up
frameRate = 3; %frames per second is hard coded unless you change this setting on the camera during acquisition
binSec = 10; % ten seconds is hard coded as the bin size unless you want to change it
totalFrames = 10800; %assumes a 1 hr long video
nanLimit = 3; % how many nan values for speed are allowed in bin before that bin is excluded
setPoint = 150; % MAY CHANGE - this is the denominator for the equation y = x/setPoint that separates roaming and dwelling based on speed vs angular speed
speedPoint = 0.05; % MAY CHANGE - speed threshold for roaming vs dwelling (above is roaming)
curvePoint = 200; % MAY CHANGE - this is the denominator for the equation y = x/curvePoint that separates roaming and dwelling based on speed vs curvature

%% Calulated from the hard coded variables
binSize = frameRate * binSec; %calculate how many frames are in ten seconds
nBins = totalFrames/binSize;
videoSec = totalFrames/frameRate;
nVideos = size(unlinkedData, 2);

%% Bin the data into 10 second intervals

% Set up a structure to store our information
rdStructure = struct('Name', cell(1,size(unlinkedData, 2)));

% Fill in the video names from the unlinkedData structure
for v = 1:nVideos
    rdStructure(v).Name= unlinkedData(v).Name;
end 

% Set up a record of the first and last frame and second in each bin
for v = 1:nVideos
    rdStructure(v).firstFrames = [1:binSize:totalFrames];
    rdStructure(v).lastFrames = [0:binSize:totalFrames];
    rdStructure(v).firstSecond = [(1/3):binSec:videoSec];
    rdStructure(v).lastSecond = [0:binSec:videoSec];
end

%% Calculate the speed in each bin

for v = 1:nVideos % for each video
    allSpeeds = []; % set up to store all the binned speed data for a single video
    for w = 1:size(unlinkedData(v).Speed, 1) % for each worm (row)
        speedData = [] ; % set up an empty vector
        speedData = unlinkedData(v).Speed(w,:); % put the worm's speed in that vector
        
        % check/fix any short or long video errors
        if size(speedData, 2) == totalFrames % if we have the right number of frames
            % do nothing
        elseif size(speedData, 2) < totalFrames % if we don't have enough frames
            speedData(end+1:totalFrames) = NaN; % pad with NaNs
        elseif size(speedData, 2) > totalFrames % if we have too many frames
            speedData = speedData(1, totalFrames);
        end
        
        % section the speedData into rows where each row has a the speed we want to bin together
        binnedSpeed = []; % start with an empty vector
        binnedSpeed = reshape(speedData, [binSize, nBins]);
        binnedSpeed = binnedSpeed'; % have to store the data this way then transform because of how the reshape function works
        rdStructure(v).binnedSpeed{w} = binnedSpeed; % store the bins in our structure
        
        % calculate the average speed in each bin
        avgSpeed = nan(1, nBins); % start with an empty vector
        for b = 1: nBins % for each bin
            if sum(isnan(binnedSpeed(b,:))) > nanLimit % don't calculate speed if you're over the limit
                avgSpeed(b)= NaN;
            else
                avgSpeed(b) = mean(binnedSpeed(b,:),'omitnan'); % otherwise take the average ignoring NaNs
            end 
        end % this should result in an avgSpeed vector where every number is the speed in 10 second interval
        
        % Save the average speeds
        rdStructure(v).avgSpeed{w} = avgSpeed; % store the bins in our structure
        
        % Add it into all the speed data
        allSpeeds = [allSpeeds; avgSpeed];
    end
    % Save the allSpeeds Data
    rdStructure(v).allSpeeds = allSpeeds;
end

%% Calculate the angular speed in each bin

for v = 1:nVideos % for each video
    allAngSpeeds = []; % set up to store all the binned speed data for a single video
    for w = 1:size(unlinkedData(v).AngSpeed, 1) % for each worm (row)
        angspeedData = [] ; % set up an empty vector
        angspeedData = unlinkedData(v).AngSpeed(w,:); % put the worm's speed in that vector
        
        % check/fix any short or long video errors
        if size(angspeedData, 2) == totalFrames % if we have the right number of frames
            % do nothing
        elseif size(angspeedData, 2) < totalFrames % if we don't have enough frames
            angspeedData(end+1:totalFrames) = NaN; % pad with NaNs
        elseif size(angspeedData, 2) > totalFrames % if we have too many frames
            angspeedData = angspeedData(1, totalFrames);
        end
        
        % section the speedData into rows where each row has a the speed we want to bin together
        binnedAngSpeed = []; % start with an empty vector
        binnedAngSpeed = reshape(angspeedData, [binSize, nBins]);
        binnedAngSpeed = binnedAngSpeed'; % have to store the data this way then transform because of how the reshape function works
        rdStructure(v).binnedAngSpeed{w} = binnedAngSpeed; % store the bins in our structure
        
        % calculate the average speed in each bin
        avgAngSpeed = nan(1, nBins); % start with an empty vector
        for b = 1: nBins % for each bin
            if sum(isnan(binnedAngSpeed(b,:))) > nanLimit % don't calculate speed if you're over the limit
                avgAngSpeed(b)= NaN;
            else
                avgAngSpeed(b) = mean(binnedAngSpeed(b,:),'omitnan'); % otherwise take the average ignoring NaNs
            end 
        end % this should result in an avgSpeed vector where every number is the speed in 10 second interval
        
        % Save the average speeds
        rdStructure(v).avgAngSpeed{w} = avgAngSpeed; % store the bins in our structure
        
        % Add it into all the speed data
        allAngSpeeds = [allAngSpeeds; avgAngSpeed];
    end
    % Save the allSpeeds Data
    rdStructure(v).allAngSpeeds = allAngSpeeds;
end

%% Calculate the curvature in each bin

for v = 1:nVideos % for each video
    allCurves = []; % set up to store all the binned speed data for a single video
    for w = 1:size(unlinkedData(v).Curvature, 1) % for each worm (row)
        curveData = [] ; % set up an empty vector
        curveData = unlinkedData(v).Curvature(w,:); % put the worm's speed in that vector
        
        % check/fix any short or long video errors
        if size(curveData, 2) == totalFrames % if we have the right number of frames
            % do nothing
        elseif size(curveData, 2) < totalFrames % if we don't have enough frames
            curveData(end+1:totalFrames) = NaN; % pad with NaNs
        elseif size(curveData, 2) > totalFrames % if we have too many frames
            curveData = curveData(1, totalFrames);
        end
        
        % section the speedData into rows where each row has a the speed we want to bin together
        binnedCurvature = []; % start with an empty vector
        binnedCurvature = reshape(curveData, [binSize, nBins]);
        binnedCurvature = binnedCurvature'; % have to store the data this way then transform because of how the reshape function works
        rdStructure(v).binnedCurvature{w} = binnedCurvature; % store the bins in our structure
        
        % calculate the average speed in each bin
        avgCurvature = nan(1, nBins); % start with an empty vector
        for b = 1: nBins % for each bin
            if sum(isnan(binnedCurvature(b,:))) > nanLimit % don't calculate speed if you're over the limit
                avgCurvature(b)= NaN;
            else
                avgCurvature(b) = mean(binnedCurvature(b,:),'omitnan'); % otherwise take the average ignoring NaNs
            end 
        end % this should result in an avgSpeed vector where every number is the speed in 10 second interval
        
        % Save the average speeds
        rdStructure(v).avgCurvature{w} = avgCurvature; % store the bins in our structure
        
        % Add it into all the speed data
        allCurves = [allCurves; avgCurvature];
    end
    % Save the allSpeeds Data
    rdStructure(v).allCurves = allCurves;
end


%% Prepare Data for Plotting

% Find out how many rows we're going to have
totalRows = 0;
for v = 1:nVideos
    addMe = size(rdStructure(v).allSpeeds, 1);
    totalRows = totalRows + addMe;
end 

% Create a matrix that has the time of each bin (reported as the middle
% time in the bin)
masterTime = repmat([(binSec/2):binSec:videoSec], totalRows, 1); % set up a matrix with the time for each bin

% Create a matrix that has all the Speed data
masterSpeed =[];
for v = 1:nVideos
    masterSpeed = cat(1, masterSpeed, rdStructure(v).allSpeeds);
end

%Create a matrix that has all the Angular Speed Data
masterAngSpeedRAW =[];
for v = 1:nVideos
    masterAngSpeedRAW = cat(1, masterAngSpeedRAW, rdStructure(v).allAngSpeeds);
end
masterAngSpeed = abs(masterAngSpeedRAW);

% Create a matrix that has all the Curvature Data
masterCurvatureRAW =[];
for v = 1:nVideos
    masterCurvatureRAW = cat(1, masterCurvatureRAW, rdStructure(v).allCurves);
end
masterCurvature = abs(masterCurvatureRAW);

%% PLOT Speed vs Time
% What to plot
X = masterTime;
Y = masterSpeed;

% Filter out NaN values
validIndices = ~isnan(X) & ~isnan(Y);
X = X(validIndices);
Y = Y(validIndices);

% Flatten the matrices to create X and Y values for the scatter plot
X = X(:);
Y = Y(:);

% Construct Scatter
 if speed_time_scatter == 1
     
     figure
     
     % Create the scatter plot
     scatter(X, Y);
     
     % Set the x-axis range
     xlim([0, 3600]);
     
     % Add labels and title
     xlabel('Time (seconds)');
     ylabel('Speed (mm/sec)');
     plotTitle = strcat(conditionName, ': speed over time');
     title(plotTitle);
     
     % Display the plot
     grid on;
     
     % Save Figure 
     figureName = strcat(conditionName, '.Speed-Time-Scatter');
     savefig(strcat(figureName,'.fig'))  % Saves as MATLAB .fig file
     saveas(gcf, strcat(figureName,'.jpg')) % Saves as a .jpg file
 end 

% Construct Heatmap
 if speed_time_heatmap == 1
     figure
     
     % Round to the nearest hundreth
     X = round(X, 2);
     Y = round(Y, 2);
     
     % Define the edges for the bins
     edgesX = unique(X);
     edgesY = unique(Y);
     
     % Ensure edges are monotonically increasing
     edgesX = sort(edgesX);
     edgesY = sort(edgesY);
     
     % Create a 2D histogram (heatmap) with the counts of each combination of X and Y values
     heatmap = histcounts2(X, Y, [edgesX; max(edgesX)+1], [edgesY; max(edgesY)+1]);
     
     % Normalize the heatmap to represent percentages
     heatmap_percentage = (heatmap / sum(heatmap(:))) * 100;
     
     % Plot the heatmap
     imagesc(edgesX, edgesY, heatmap_percentage');
     colorbar;
     colormap('hot');
     caxis([0, max(heatmap_percentage(:))]);
     
     % Set YDir to 'normal' to have the lowest Y value at the bottom
     set(gca, 'YDir', 'normal');
     
     % Set the x-axis range
     xlim([0, 3600]);
     
     % Add labels and title
     xlabel('Time (seconds)');
     ylabel('Speed (mm/sec)');
     plotTitle = strcat(conditionName, ': speed over time');
     title(plotTitle);
     
     % Display the plot
     grid on;
     
     % Plot dividing slope
     hold on;
     x = 0:3600;  % Match your x-axis range
     y = ones(size(x)) * speedPoint;  % Constant threshold line at speedPoint
     plot(x, y, 'r-', 'LineWidth', 2);
     hold off;
     
     % Save Figure 
     figureName = strcat(conditionName, '.Speed-Time-Heatmap');
     savefig(strcat(figureName,'.fig'))  % Saves as MATLAB .fig file
     saveas(gcf, strcat(figureName,'.jpg')) % Saves as a .jpg file
     
 end 

%% PLOT Angular Speed vs Time
% What to plot
X = masterTime;
Y = masterAngSpeed;

% Filter out NaN values
validIndices = ~isnan(X) & ~isnan(Y);
X = X(validIndices);
Y = Y(validIndices);
% Construct Scatter

figure

% Create the scatter plot
if ang_time_scatter == 1
    scatter(X, Y);
    
    
    % Set the x-axis and y-axis range
    xlim([0, 3600]);
    %ylim([0,180]);
    
    % Add labels and title
    xlabel('Time (seconds)');
    ylabel('Angular Speed (deg/sec)');
    plotTitle = strcat(conditionName, ': angular speed over time');
    title(plotTitle);
    
    % Display the plot
    grid on;
    
    % Save Figure
    figureName = strcat(conditionName, '.Ang-Time-Scatter');
    savefig(strcat(figureName,'.fig'))  % Saves as MATLAB .fig file
    saveas(gcf, strcat(figureName,'.jpg')) % Saves as a .jpg file
end 

% Construct Heatmap
 if ang_time_heatmap == 1
     figure
     
     % Filter out NaN values
     validIndices = ~isnan(X) & ~isnan(Y);
     X = X(validIndices);
     Y = Y(validIndices);
     
     % Round to the nearest hundreth
     X = round(X, 2);
     Y = round(Y, 2);
     
     % Define the edges for the bins
     edgesX = unique(X);
     edgesY = unique(Y);
     
     % Ensure edges are monotonically increasing
     edgesX = sort(edgesX);
     edgesY = sort(edgesY);
     
     % Create a 2D histogram (heatmap) with the counts of each combination of X and Y values
     heatmap = histcounts2(X, Y, [edgesX; max(edgesX)+1], [edgesY; max(edgesY)+1]);
     
     % Normalize the heatmap to represent percentages
     heatmap_percentage = (heatmap / sum(heatmap(:))) * 100;
     
     % Plot the heatmap
     imagesc(edgesX, edgesY, heatmap_percentage');
     colorbar;
     colormap('hot');
     caxis([0, max(heatmap_percentage(:))]);
     
     % Set YDir to 'normal' to have the lowest Y value at the bottom
     set(gca, 'YDir', 'normal');
     
     % Set the x-axis and y-axis range
     xlim([0, 3600]);
     %ylim([0,180]);
     
     % Add labels and title
     xlabel('Time (seconds)');
     ylabel('Angular Speed (deg/sec)');
     plotTitle = strcat(conditionName, ': angular speed over time');
     title(plotTitle);
     
     % Display the plot
     grid on;
     
     % Save Figure 
     figureName = strcat(conditionName, '.Ang-Time-Scatter');
     savefig(strcat(figureName,'.fig'))  % Saves as MATLAB .fig file
     saveas(gcf, strcat(figureName,'.jpg')) % Saves as a .jpg file
 end 
%% PLOT Speed vs Angular Speed
% What to plot
X = masterAngSpeed;
Y = masterSpeed;

% Filter out NaN values
validIndices = ~isnan(X) & ~isnan(Y);
X = X(validIndices);
Y = Y(validIndices);

% Construct Scatter
if speed_ang_scatter == 1
    figure
    
    % Create the scatter plot
    scatter(X, Y);
    
    % Set the x-axis range
    %xlim([0, 180]);
    
    % Add labels and title
    xlabel('Angular Speed (deg/sec)');
    ylabel('Speed (mm/sec)');
    plotTitle = strcat(conditionName, ': angular speed vs speed');
    title(plotTitle);
    
    % Display the plot
    grid on;
    
    % Save Figure
    figureName = strcat(conditionName, '.Speed-Ang-Scatter');
    savefig(strcat(figureName,'.fig'))  % Saves as MATLAB .fig file
    saveas(gcf, strcat(figureName,'.jpg')) % Saves as a .jpg file
end 

% Construct Heatmap
if speed_ang_heatmap == 1
    figure
    
    % Filter out NaN values
    validIndices = ~isnan(X) & ~isnan(Y);
    X = X(validIndices);
    Y = Y(validIndices);
    
    % Round to the nearest hundreth
    X = round(X, 2);
    Y = round(Y, 2);
    
    % Define the edges for the bins
    edgesX = unique(X);
    edgesY = unique(Y);
    
    % Ensure edges are monotonically increasing
    edgesX = sort(edgesX);
    edgesY = sort(edgesY);
    
    % Create a 2D histogram (heatmap) with the counts of each combination of X and Y values
    heatmap = histcounts2(X, Y, [edgesX; max(edgesX)+1], [edgesY; max(edgesY)+1]);
    
    % Normalize the heatmap to represent percentages
    heatmap_percentage = (heatmap / sum(heatmap(:))) * 100;
    
    % Plot the heatmap
    imagesc(edgesX, edgesY, heatmap_percentage');
    colorbar;
    colormap('hot');
    caxis([0, max(heatmap_percentage(:))]);
    
    % Set YDir to 'normal' to have the lowest Y value at the bottom
    set(gca, 'YDir', 'normal');
    
    % Set the x-axis range
    %xlim([0, 180]);
    
    % Add labels and title
    xlabel('Angular Speed (deg/sec)');
    ylabel('Speed (mm/sec)');
    plotTitle = strcat(conditionName, ': angular speed over time');
    title(plotTitle);
    
    % Display the plot
    grid on;
    
    % Plot dividing slope
    hold on;
    x = 0:180;  % Assuming your x-axis goes to 180 deg/sec
    y = x/setPoint;  %(default 1/450, range 1/270 - 1/550 in Flavell 2013)
    plot(x, y, 'r-', 'LineWidth', 2);
    hold off;
    
    % Label the heatmap key
    c = colorbar;
    c.Label.String = 'Fraction of data points';
    c.Label.FontSize = 12;
    c.Label.Rotation = 270;
    
    % Adjust the axes position to make room for colorbar
    c.Label.Position = c.Label.Position + [1 0 0];  % Move label slightly to the right
    
    % Save Figure
    figureName = strcat(conditionName, '.Speed-Ang-Heatmap');
    savefig(strcat(figureName,'.fig'))  % Saves as MATLAB .fig file
    saveas(gcf, strcat(figureName,'.jpg')) % Saves as a .jpg file
end 

%% PLOT Speed vs Curvature
% What to plot
X = masterCurvature;
Y = masterSpeed;

% Filter out NaN values
validIndices = ~isnan(X) & ~isnan(Y);
X = X(validIndices);
Y = Y(validIndices);

% Construct Scatter
if speed_curve_scatter == 1
    figure
    
    % Create the scatter plot
    scatter(X, Y);
    
    % Set the x-axis range
    %xlim([0, 180]);
    
    % Add labels and title
    xlabel('Curvature(deg)');
    ylabel('Speed (mm/sec)');
    plotTitle = strcat(conditionName, ': speed vs curvature');
    title(plotTitle);
    
    % Display the plot
    grid on;
    
    % Save Figure
    figureName = strcat(conditionName, '.Speed-Curve-Scatter');
    savefig(strcat(figureName,'.fig'))  % Saves as MATLAB .fig file
    saveas(gcf, strcat(figureName,'.jpg')) % Saves as a .jpg file
end 

% Construct Heatmap
if speed_curve_heatmap == 1
    figure
    
    % Filter out NaN values
    validIndices = ~isnan(X) & ~isnan(Y);
    X = X(validIndices);
    Y = Y(validIndices);
    
    % Round to the nearest hundreth
    X = round(X, 2);
    Y = round(Y, 2);
    
    % Define the edges for the bins
    edgesX = unique(X);
    edgesY = unique(Y);
    
    % Ensure edges are monotonically increasing
    edgesX = sort(edgesX);
    edgesY = sort(edgesY);
    
    % Create a 2D histogram (heatmap) with the counts of each combination of X and Y values
    heatmap = histcounts2(X, Y, [edgesX; max(edgesX)+1], [edgesY; max(edgesY)+1]);
    
    % Normalize the heatmap to represent percentages
    heatmap_percentage = (heatmap / sum(heatmap(:))) * 100;
    
    % Plot the heatmap
    imagesc(edgesX, edgesY, heatmap_percentage');
    colorbar;
    colormap('hot');
    caxis([0, max(heatmap_percentage(:))]);
    
    % Set YDir to 'normal' to have the lowest Y value at the bottom
    set(gca, 'YDir', 'normal');
    
    % Set the x-axis range
    %xlim([0, 180]);
    
    % Add labels and title
    xlabel('Curvature (deg)');
    ylabel('Speed (mm/sec)');
    plotTitle = strcat(conditionName, ': speed vs curvature');
    title(plotTitle);
    
    % Display the plot
    grid on;

    % Plot dividing slope
    hold on;
    x = 0:180;  % Assuming your x-axis goes to 180 deg/sec
    y = x/curvePoint;  %(default 1/200 extrapolated from Ben Arous 2009)
    plot(x, y, 'r-', 'LineWidth', 2);
    hold off;
    
        
    % Save Figure
    figureName = strcat(conditionName, '.Speed-Curve-Heatmap');
    savefig(strcat(figureName,'.fig'))  % Saves as MATLAB .fig file
    saveas(gcf, strcat(figureName,'.jpg')) % Saves as a .jpg file
end 


%% Plot the speed over time by plate

plateSpeeds = nan(nVideos, nBins);
for v = 1:nVideos
    plateSpeeds(v,:)= mean(rdStructure(v).allSpeeds,'omitnan');
end

if speed_time_plate == 1
    figure
    
    plot(masterTime(1,:), plateSpeeds)
    
    
    % Set the x-axis range
    xlim([0, 3600]);
    
    % Add labels and title
    xlabel('Time (sec)');
    ylabel('Speed (mm/sec)');
    plotTitle = strcat(conditionName, ': speed by plate');
    title(plotTitle);
    
    % Display the plot
    grid on;
    
    % Save Figure
    figureName = strcat(conditionName, '.Speed-Time-Plate');
    savefig(strcat(figureName,'.fig'))  % Saves as MATLAB .fig file
    saveas(gcf, strcat(figureName,'.jpg')) % Saves as a .jpg file
end 

%% Calculate Roaming and Dwelling Based on Angular Speed and Velocity

% Reshape data to vectors, removing NaNs
validIdx = ~isnan(masterSpeed) & ~isnan(masterAngSpeed);
speed = masterSpeed(validIdx);
angSpeed = masterAngSpeed(validIdx);

% Allow slope adjustment (default 1/450, range 1/270 - 1/550)
slope = 1/setPoint;

% Classify into behaviors
observations = zeros(size(speed));
for i = 1:length(speed)
    if speed(i) > angSpeed(i)/setPoint  % Points above line are roaming
        observations(i) = 1; % Roam
    else
        observations(i) = 2; % Dwell
    end
end

% Define transition probabilities
tr = [0.9412, 0.0588;   
     0.0168, 0.9832];  

% Define emission probabilities
em = [0.9037, 0.0963;   
     0.0196, 0.9804];  

% Get state sequence
estimatedStates = hmmviterbi(observations, tr, em);

% Reshape back to original matrix
stateMatrix = nan(size(masterSpeed));
stateMatrix(validIdx) = estimatedStates;

% Calculate fractions
fractionDwelling = sum(estimatedStates == 2) / length(estimatedStates);
fractionRoaming = sum(estimatedStates == 1) / length(estimatedStates);

% Create stacked bar plot
if roam_dwell_speed_ang == 1
    figure;
    bar(1:1, [fractionDwelling fractionRoaming], 0.5, 'stack');
    colormap([0 0 0.8; 1 0.4 0]);
    ylim([0 1]);
    set(gca, 'XTick', 1, 'XTickLabel', 'WT');
    ylabel('fraction of time spent in state');
    legend({'Dwelling', 'Roaming'}, 'Location', 'eastoutside');
    plotTitle = strcat(conditionName, ': Roaming vs. Dwelling');
    title(plotTitle);
    
    % Save Figure
    figureName = strcat(conditionName, '.RD-Speed_Ang');
    savefig(strcat(figureName,'.fig'))  % Saves as MATLAB .fig file
    saveas(gcf, strcat(figureName,'.jpg')) % Saves as a .jpg file
end 

rdmatrixAngSpeed = stateMatrix;
stateStruct = struct;
stateStruct(1).Name = 'Ang Speed v Velocity';
stateStruct(1).Dwelling = fractionDwelling;
stateStruct(1).Roaming = fractionRoaming;

%% Calulate Roaming and Dwelling based only on Speed

% Use just speed data
validIdx = ~isnan(masterSpeed);
speed = masterSpeed(validIdx);

% Classify based on speed threshold
observations = zeros(size(speed));
observations(speed > speedPoint) = 1; % Roam
observations(speed <= speedPoint) = 2; % Dwell

% HMM
estimatedStates = hmmviterbi(observations, tr, em);

% Calculate fractions
fractionDwelling = sum(estimatedStates == 2) / length(estimatedStates);
fractionRoaming = sum(estimatedStates == 1) / length(estimatedStates);

% Plot
if roam_dwell_only_speed == 1
    figure;
    bar(1:1, [fractionDwelling fractionRoaming], 0.5, 'stack');
    colormap([0 0 0.8; 1 0.4 0]);
    ylim([0 1]);
    set(gca, 'XTick', 1, 'XTickLabel', 'WT');
    ylabel('fraction of time spent in state');
    legend({'Dwelling', 'Roaming'}, 'Location', 'eastoutside');
    plotTitle = strcat(conditionName, ': Roaming vs. Dwelling');
    title(plotTitle);
    
    % Save Figure 
    figureName = strcat(conditionName, '.RD-Speed_Ang');
    savefig(strcat(figureName,'.fig'))  % Saves as MATLAB .fig file
    saveas(gcf, strcat(figureName,'.jpg')) % Saves as a .jpg file
end

rdmatrixSpeedOnly = stateMatrix;
stateStruct(2).Name = 'Speed Only';
stateStruct(2).Dwelling = fractionDwelling;
stateStruct(2).Roaming = fractionRoaming;

%% Calculate Roaming and Dwelling Based on Speed vs Curvature 

% Reshape data to vectors, removing NaNs
validIdx = ~isnan(masterSpeed) & ~isnan(masterCurvature);
speed = masterSpeed(validIdx);
curvature = masterCurvature(validIdx);

% Allow slope adjustment (default 1/450, range 1/270 - 1/550)
slope = 1/curvePoint;

% Classify into behaviors
observations = zeros(size(speed));
for i = 1:length(speed)
    if speed(i) > curvature(i)/curvePoint  % Points above line are roaming
        observations(i) = 1; % Roam
    else
        observations(i) = 2; % Dwell
    end
end

% Define transition probabilities
tr = [0.9412, 0.0588;   
     0.0168, 0.9832];  

% Define emission probabilities
em = [0.9037, 0.0963;   
     0.0196, 0.9804];  

% Get state sequence
estimatedStates = hmmviterbi(observations, tr, em);

% Reshape back to original matrix
stateMatrix = nan(size(masterSpeed));
stateMatrix(validIdx) = estimatedStates;

% Calculate fractions
fractionDwelling = sum(estimatedStates == 2) / length(estimatedStates);
fractionRoaming = sum(estimatedStates == 1) / length(estimatedStates);

% Create stacked bar plot
if roam_dwell_speed_curve == 1
    figure;
    bar(1:1, [fractionDwelling fractionRoaming], 0.5, 'stack');
    colormap([0 0 0.8; 1 0.4 0]);
    ylim([0 1]);
    set(gca, 'XTick', 1, 'XTickLabel', 'WT');
    ylabel('fraction of time spent in state');
    legend({'Dwelling', 'Roaming'}, 'Location', 'eastoutside');
    plotTitle = strcat(conditionName, ': Roaming vs. Dwelling');
    title(plotTitle);
    
    % Save Figure
    figureName = strcat(conditionName, '.RD-Speed_Curve');
    savefig(strcat(figureName,'.fig'))  % Saves as MATLAB .fig file
    saveas(gcf, strcat(figureName,'.jpg')) % Saves as a .jpg file
end

rdmatrixSpeedCurve = stateMatrix;
stateStruct(3).Name = 'Ang Speed v Velocity';
stateStruct(3).Dwelling = fractionDwelling;
stateStruct(3).Roaming = fractionRoaming;

%% Close Figures
if closeThem == 1
    close all
end 
    
%% Store and Save Variables

% Save rdStructure
varName = strcat(conditionName, '.binnedData.mat');
save(varName', 'rdStructure')

% Store and save processed data
processedData = struct;
processedData.masterTime = masterTime;
processedData.masterSpeed = masterSpeed;
processedData.masterAngSpeedRAW = masterAngSpeedRAW;
processedData.masterAngSpeed = masterAngSpeed;
processedData.masterCurvature = masterCurvature;
processedData.plateSpeeds = plateSpeeds;
processedData.setPoint = setPoint;
processedData.speedPoint = speedPoint;
processedData.curvePoint = curvePoint;
processedData.AngSpeedSpeedRD = rdmatrixAngSpeed;
processedData.SpeedOnlyRD = rdmatrixSpeedOnly;
processedData.SpeedCurveRD = rdmatrixSpeedCurve;
processedData.Fractions = stateStruct;

varName = strcat(conditionName, '.processedData.mat');
save(varName', 'processedData')

