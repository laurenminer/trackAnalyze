function mat = get_mat_from_tracks(tracks,property)
%{
Inputs
- tracks: linkedTracks structure (or chxTracks) loaded from the .mat file
- property: string corresponding to a field of tracks

Output
- mat: time-aligned matrix of property for all tracks
    Matrix of size length(linkedTracks) x video length in frames
    Each row i is the value of tracks(i).(property) at that video frame
    NaN everywhere that there is no data.
%}

tracks = fix_tracks_allfields(tracks); % don't return interpolated data
n_tracks = length(tracks);

% Video length in frames is not stored directly so find the largest
% frame value in Frames field.
n_frames = 0;
for i = 1:n_tracks
    last_frame = tracks(i).Frames(end);
    if last_frame > n_frames
        n_frames = last_frame;
    end
end

% Construct time-aligned matrix
mat = NaN(n_tracks,n_frames); % All "missing" frames are NaNs
for j = 1:length(tracks)
    nonzero_frames = tracks(j).Frames; % Get frame indices in which there is speed data
    mat(j,nonzero_frames) = tracks(j).(property); % Map the speed data onto those indices (rest remain NaNs)
end

return