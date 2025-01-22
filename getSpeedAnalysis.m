function getSpeedAnalysis(unlinkedDataName)


% Calculates the speed in 10 second intervals
% Intervals are pre-set as consecutive blocks of 30 frames
% Calculation is not performed if there are more than 3 frames missing

%% Load Data
load(unlinkedDataName)

%% Variables/Housekeeping

% Hard coded variables determined by experimental set up
frameRate = 3; %frames per second is hard coded unless you change this setting on the camera during acquisition
binSec = 10; % ten seconds is hard coded as the bin size unless you want to change it
totalFrames = 10800; %assumes a 1 hr long video
nanLimit = 3; % how many nan values for speed are allowed in bin before that bin is excluded
setPoint = 200; % MAY CHANGE - this is the denominator for the equation y = x/setPoint that separates roaming and dwelling
speedPoint = 0.05; % MAY CHANGE - speed threshold for roaming vs dwelling

%Calulated from the hard coded variables
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
    for w = 1:size(unlinkedData(v).Speed, 1) % for each worm (row)
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

%Create a matriz that has all the Angular Speed Data
masterAngSpeedRAW =[];
for v = 1:nVideos
    masterAngSpeedRAW = cat(1, masterAngSpeedRAW, rdStructure(v).allAngSpeeds);
end
masterAngSpeed = abs(masterAngSpeedRAW);

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

% % Construct Scatter
% 
% figure
% 
% % Create the scatter plot
% scatter(X, Y);
% 
% % Set the x-axis range
% xlim([0, 3600]);
% 
% % Add labels and title
% xlabel('Time (seconds)');
% ylabel('Speed (mm/sec)');
% title('Worm speed over time');
% 
% % Display the plot
% grid on;

% Construct Heatmap

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
title('Worm speed over time');

% Display the plot
grid on;

% Plot dividing slope
hold on;
x = 0:3600;  % Match your x-axis range
y = ones(size(x)) * speedPoint;  % Constant threshold line at speedPoint
plot(x, y, 'r-', 'LineWidth', 2);
hold off;

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
scatter(X, Y);


% Set the x-axis and y-axis range
xlim([0, 3600]);
%ylim([0,180]);

% Add labels and title
xlabel('Time (seconds)');
ylabel('Angular Speed (deg/sec)');
title('Angular speed over time');

% Display the plot
grid on;

% Construct Heatmap

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
title('Angular speed over time');

% Display the plot
grid on;

%% PLOT Speed vs Angular Speed
% What to plot
X = masterAngSpeed;
Y = masterSpeed;

% Filter out NaN values
validIndices = ~isnan(X) & ~isnan(Y);
X = X(validIndices);
Y = Y(validIndices);

% Construct Scatter

figure

% Create the scatter plot
scatter(X, Y);

% Set the x-axis range
%xlim([0, 180]);

% Add labels and title
xlabel('Angular Speed (deg/sec)');
ylabel('Speed (mm/sec)');
title('Worm speed over time');

% Display the plot
grid on;

% Construct Heatmap

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
title('Speed and Angular Speed');

% Display the plot
grid on;

% Plot dividing slope
hold on;
x = 0:180;  % Assuming your x-axis goes to 180 deg/sec
y = x/setPoint;  %(default 1/450, range 1/270 - 1/550 in Flavell 2013)
plot(x, y, 'r-', 'LineWidth', 2);
hold off;

%% Plot the speed over time by plate

plateSpeeds = nan(nVideos, nBins);
for v = 1:nVideos
    plateSpeeds(v,:)= mean(rdStructure(v).allSpeeds,'omitnan');
end

figure 

plot(masterTime(1,:), plateSpeeds)


% Set the x-axis range
xlim([0, 3600]);

% Add labels and title
xlabel('Time (sec)');
ylabel('Speed (mm/sec)');
title('Speed over time');

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
figure;
bar(1:1, [fractionDwelling fractionRoaming], 0.5, 'stack');
colormap([0 0 0.8; 1 0.4 0]);
ylim([0 1]);
set(gca, 'XTick', 1, 'XTickLabel', 'WT');
ylabel('fraction of time spent in state');
legend({'Dwelling', 'Roaming'}, 'Location', 'eastoutside');

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
figure;
bar(1:1, [fractionDwelling fractionRoaming], 0.5, 'stack');
colormap([0 0 0.8; 1 0.4 0]);
ylim([0 1]);
set(gca, 'XTick', 1, 'XTickLabel', 'WT');
ylabel('fraction of time spent in state');
legend({'Dwelling', 'Roaming'}, 'Location', 'eastoutside');