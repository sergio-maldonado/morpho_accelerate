%accelerate_simulation
%%% FUNCTION TO ACCELERATE XBEACH SIMULATION %%% 
% Author: Ellie Newell, University of Southampton

function accelerate_simulation()
    %% Set up
    % Declare global variables
    global x
    global z_initial
    global tintg
    global simulation_time
    global extrapolate_time
    global sample_interval
    global extrapolation_fit
    global run_XB
    global new_tstart
    
    % Run XBeach simulation if required
    if run_XB == 1 
        system('xbeach.exe'); % run XBeach simulation 
    end
    %% Extrapolate x points
    % Set up matrices of sample points
    time_matrix = [1:tintg:(simulation_time + 1 - new_tstart)]; % time matrix of simulation
    sample_matrix = [1:sample_interval:length(x)]; % matrix of sample x points
    sample_x_points = zeros(1,length(sample_matrix)); % values of sample x points
    sample_z_points = zeros(1,length(sample_matrix)); % values of sample z points
    for xpoint = 1:length(sample_matrix) % fill sample matrices
        sample_x_points(xpoint) = x(sample_matrix(xpoint));
        sample_z_points(xpoint) = z_initial(sample_matrix(xpoint));
    end
    if sample_x_points(length(sample_matrix)) ~= x(length(x)) % add final point to sample matrix (if length(x) not multiple of spacing interval)
        sample_x_points(length(sample_matrix)+1) = x(length(x));
        sample_z_points(length(sample_matrix)+1) = z_initial(length(x));
        sample_matrix(1, length(sample_matrix)+1) = length(x); % add final xpoint to sample matrix
    end
    
    % Bed change at each time for sample points
    zb = xb_read_dat('zb.dat'); % read .dat output into a MATLAB structure
    t = length(zb.data(2).value(:,:,1)); % find the simulation time t
    bed_diff = zeros(length(sample_x_points), t); % create empty matrix to populate with bed differences
    for time = 1:t % extract data from zb structure output for each time 
        zb_data = zb.data(2).value(time,:,:); % extract data for time i
        for xpoint = 1:length(sample_x_points) % find difference between initial bed and bed depth at time i 
            bed_diff(xpoint,time) = -1*(sample_z_points(xpoint) - zb_data(sample_matrix(xpoint))); % write bed difference to matrix 
        end
    end
    
    % Extrapolate sample points based on bed change curve 
    sample_predicted_change = zeros(1, length(x)); % empty matrix for predicted changes
    sample_extrapolated_profile = zeros(1,length(sample_matrix)); % empty matrix for new profile
    for xpoint = 1:length(sample_x_points) % for each point
        sample_extrapolate_trend = polyval(polyfit(time_matrix, bed_diff(xpoint,time_matrix),extrapolation_fit), simulation_time + extrapolate_time); % evaluate line fit to bed change over time 
        sample_predicted_change(1, xpoint) = sample_extrapolate_trend; % add extrapolate_trend to predicted change matrix
        sample_extrapolated_profile(1, xpoint) = sample_z_points(xpoint) + sample_predicted_change(xpoint); % calculate new bed position
    end
    
    % Fix boundary points
    sample_extrapolated_profile(1,1) = z_initial(1); % offshore boundary
    sample_extrapolated_profile(1,length(sample_x_points)) = z_initial(length(x)); % onshore boundary
    
    % Interpolate with spline between sample points 
    global extrapolated_profile
    extrapolated_profile = interp1(sample_x_points,sample_extrapolated_profile,x,'spline');
    
    %% Update file
    % Update bed.dep file
    depfile = fopen('bed.dep', 'w'); % open bed.dep to write
    fprintf(depfile, '%d ', extrapolated_profile); % overwrite with predicted_change
    fclose(depfile); % close bed.dep
end
