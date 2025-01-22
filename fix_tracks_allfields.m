function fixedTracks = fix_tracks_allfields(linkedTracks)
%{
Takes in a linkedTracks structure and removes data that is just made up/
interpolated during time points where the tracker lost the worm. This
happens in many fields and can mess up measures like speed and position.
%}

fixedTracks = linkedTracks;

fieldNames = fieldnames(linkedTracks);

% Loop through fields
for fieldIdx = 1:length(fieldNames)
    field = fieldNames{fieldIdx};
    
    % The frame and time corresponding to the entire track are still useful
    % and important for indexing, even if the data recorded during that frame is deleted
    if strcmp(field,'Frames') || strcmp(field,'Time')
        continue
    end

    testVar = linkedTracks(1).(field);
    % Only execute if field contains a vector of one value for each frame:
    if isa(testVar,"double") || isa(testVar,"single")
        if length(testVar) == length(linkedTracks(1).Frames)

            % If it meets these conditions, loop through each track and
            % replace interpolated data with NaN
            for track = 1:length(linkedTracks)    
                currVar = linkedTracks(track).(field);
                % midbody_angle is not interpolated and has NaNs where there's no real data
                currVar(isnan(linkedTracks(track).midbody_angle)) = NaN;
                fixedTracks(track).(field) = currVar;
            end

        end
    end
end

return