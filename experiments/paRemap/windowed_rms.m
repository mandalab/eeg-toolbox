function result = windowed_rms(data,points_per_win,points_per_slide)
    %This computes the moving window average of data which is time X
    %channels
    total_points = size(data,1);
    n_channels = size(data,2);
    points_per_overlap = points_per_win - points_per_slide;
    n_windows = floor((total_points - points_per_overlap)/points_per_slide);
    result = zeros(n_windows,n_channels);
    for i = 1:n_windows
        start = (i-1)*points_per_slide + 1;
        result(i,:) = rms_norm(data(start:start+points_per_win-1,:));
    end
end