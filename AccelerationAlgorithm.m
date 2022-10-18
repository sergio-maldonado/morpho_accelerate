%%% ACCELERATION ALGORITHM %%%
% Author: Ellie Newell, University of Southampton

%% Clear MATLAB environment and start timer

clc; clear variables; close all
tic; % start timer

%% Extract data from params.txt and jonswap.txt files

% Using extract_params.m function
D50 = extract_params('D50'); % extract D50 value
zs0 = extract_params('zs0'); % extract initial water level value
global posdwn % declare as global variable
posdwn = extract_params('posdwn'); % extract posdwn value 
global original_tstop % declare as global variable
original_tstop = extract_params('tstop'); % extract tstop value
global tintg % declare as global variable
tintg = extract_params('tintg'); % extract tintg value

% Using  extract_jons.m
H = extract_jons('Hm0'); % extract Hm0 value
fp = extract_jons('fp'); % extract fp value
T = 1/fp; % wave period = 1/frequency

%% Import beach profile data 

global x % declare as global variable
x = importdata('x.grd'); % import x coordinates from x.grd file
global z_initial % declare as global variable
z_initial = importdata('bed.dep'); % import z coordinates from bed.dep file
if posdwn == 1 % if posdwn is 1
    z_initial = -1.*z_initial; % flip profile to account for different positive directions
end

%% Calculate non-dimensional time parameter

% Find median bed depth
water_depth = zeros(1,length(x)); % empty initial water depth matrix
for xpoint = 1:length(x) % for all xpoints 
   depth = zs0 - z_initial(xpoint); % calculate water depth
   if depth > 0 % if xpoint is below initial water level
       water_depth(xpoint) = depth; % add depth value to matrix
   end 
end
median_depth = median(water_depth); % find median water depth

% Calculate value of non-dimensional parameter
nond_param = (D50*median_depth)/(H^2); % larger nond_param = longer extrapolation times

%% Calculate tstart and morstart

g = 9.81; % gravitational acceleration (ms^-2)
wavelength = (g*(T^2))/(2*pi()); % wavelength (m)
wavespeed = sqrt(((g*wavelength)/(2*pi()))*tanh(((2*pi())*median_depth)/wavelength)); % wavespeed (ms^-1)
% find x point closest to initial water level (first zero depth)
for xpoint = 1:length(x) % for each x point
    if water_depth(xpoint) ~= 0 % if the water depth is not zero
        propagation_distance = x(xpoint+1); % propagation distance = x distance of next point
    end
end

propagation_time = 2*(propagation_distance/wavespeed); % time for waves to propagate to shore and reflect (s)
global new_tstart % declare as global variable
new_tstart = ceil(propagation_time); % round propagation_time up to nearest integer

%% Determine simulation extrapolation time, extrapolation interval

% Time XBeach should be run for 
simulation_time_1 = 3000;
simulation_time_2 = 2000;

% Time to extrapolate changes over 
extrapolate_time_1 = 5000;
extrapolate_time_2 = 6000;

% Sample interval for x points in extrapolation
sample_interval_1 = 5;
sample_interval_2 = 1;

% Extrapolation curve fit (1 = linear, 2 = quadratic)
extrapolation_fit_1 = 1;
extrapolation_fit_2 = 2;

% Compare non-dimensional parameter to thresholds to find extrapolation time
global simulation_time % declare as global variable
global extrapolate_time % declare as global variable
global sample_interval % declare as global variable
global extrapolation_fit % declare as global variable
if nond_param < 1e-04 % if non-dimensional parameter less than 1x10^-4
    simulation_time = simulation_time_1; % simulation_time set to value 1
    extrapolate_time = extrapolate_time_1; % extrapolate_time set to value 1
    sample_interval = sample_interval_1; % sample interval set to value 1
    extrapolation_fit = extrapolation_fit_1; % extrapolation fit set to value 1
elseif 1e-04 <= nond_param && nond_param < 1e-03 % if non-dimensional parameter between 1x10^-4 and 1x10^-3
    simulation_time = simulation_time_1; % simulation_time set to value 1
    extrapolate_time = extrapolate_time_1; % extrapolate_time set to value 1
    sample_interval = sample_interval_2; % sample interval set to value 2
    extrapolation_fit = extrapolation_fit_1; % extrapolation fit set to value 1
elseif 1e-03 <= nond_param && nond_param < 1e-02 % if non-dimensional parameter between 1x10^-3 and 1x10^-2
    simulation_time = simulation_time_2; % simulation_time set to value 1
    extrapolate_time = extrapolate_time_1; % extrapolate_time set to value 1
    sample_interval = sample_interval_2; % sample interval set to value 3
    extrapolation_fit = extrapolation_fit_1; % extrapolation fit set to value 2
elseif nond_param >= 1e-02 % if non-dimensional parameter greater than 1x10^-2
    simulation_time = simulation_time_1; % simulation_time set to value 2
    extrapolate_time = extrapolate_time_2; % extrapolate_time set to value 2
    sample_interval = sample_interval_2; % sample interval set to value 3
    extrapolation_fit = extrapolation_fit_1; % extrapolation fit set to value 1
end

%% Overwrite values of tstart, morfac and tstop

params = fopen('params.txt', 'r+'); % open params.txt file to read and write
read_params = textscan(params, '%s', 'CommentStyle', '%', 'TextType', 'string'); % read params into a cell array

% Find and overwrite tstart
tstart_pat = contains(read_params{1,1}, 'tstart'); % array of 1 or 0 if params_data has variable
tstart_in_params = find(tstart_pat); % index of variable location
if tstart_in_params > 0 % if tstart specified in params.txt
    read_params{1,1}(tstart_in_params+2, 1) = new_tstart; % rewrite read_params
else % if not
    read_params{1,1}(length(read_params{1,1})+1,1) = 'tstart'; % add tstart to read_params
    read_params{1,1}(length(read_params{1,1})+1,1) = '='; % add = delimiter to read_params
    read_params{1,1}(length(read_params{1,1})+1,1) = new_tstart; % add tstart value to read_params
end

% Find and overwrite morstart
morstart_pat = contains(read_params{1,1}, 'morstart'); % array of 1 or 0 if params_data has variable
morstart_in_params = find(morstart_pat); % index of variable location
new_morstart = new_tstart; % set morstart = tstart
if morstart_in_params > 0 % if morstart specified in params.txt
    read_params{1,1}(morstart_in_params+2, 1) = new_morstart; % rewrite read_params
else % if not
    read_params{1,1}(length(read_params{1,1})+1,1) = 'morstart'; % add morstart to read_params
    read_params{1,1}(length(read_params{1,1})+1,1) = '='; % add = delimiter to read_params
    read_params{1,1}(length(read_params{1,1})+1,1) = new_morstart; % add morstart value to read_params
end

% Find and overwrite morfac
morfac_pat = contains(read_params{1,1}, 'morfac'); % array of 1 or 0 if params_data has variable
morfac_in_params = find(morfac_pat); % index of variable location
new_morfac = num2cell(1); % set morfac to 1 (no acceleration by XBeach)
if morfac_in_params > 0 % if morfac specified in params.txt
    read_params{1,1}(morfac_in_params+2, 1) = new_morfac; % rewrite read_params
else % if not
    read_params{1,1}(length(read_params{1,1})+1,1) = 'morfac'; % add morfac to read_params
    read_params{1,1}(length(read_params{1,1})+1,1) = '='; % add = delimiter to read_params
    read_params{1,1}(length(read_params{1,1})+1,1) = new_morfac; % add morfac value to read_params
end

% Find and replace tstop
tstop_pat = contains(read_params{1,1}, 'tstop'); % array of 1 or 0 if params_data has variable
tstop_in_params = find(tstop_pat); % index of variable location
new_tstop = num2cell(simulation_time); % new value of variable is simulation_time
read_params{1,1}(tstop_in_params+2, 1) = new_tstop; % rewrite read_params 

% Find and replace tintg
tintg_pat = contains(read_params{1,1}, 'tintg'); % array of 1 or 0 if params_data has variable
tintg_in_params = find(tintg_pat); % index of variable location
new_tintg = num2cell(1); % set tintg to 1s
read_params{1,1}(tintg_in_params+2, 1) = new_tintg; % rewrite read_params

% Add new values to params array
new_params = string.empty;
a = 1; % counter in for loop
for variable = 1:length(read_params{1,1})
    variable_name = rem(variable,3); % calculate if 'data' is divisible by 3 (if yes it is a variable name, if no it is a value or =)
    if variable_name == 1 % if 'data' is odd
        new_params(a,1) = read_params{1,1}(variable); % find variable and write to array
        new_params(a,2) = read_params{1,1}(variable+1); % equals sign and write to array
        new_params(a,3) = read_params{1,1}(variable+2); % find value and write to array
        a = a + 1; % increase count
    end
end

writematrix(new_params, 'params.txt', 'Delimiter', ' '); % overwrite params file

%% Run XBeach simulation and extrapolate bed changes

current_time = 0 % track time in simulation
XB_iteration_count = 0 % count number of times XBeach runs
global run_XB % declare global variable

while current_time < original_tstop % while time in simulation less than original tstop
    if original_tstop - current_time - simulation_time - extrapolate_time > 0 % if current time + simulation time + extrapolation time less than original tstop 
        run_XB = 1; % tell accelerate_simulation function to run XBeach
        accelerate_simulation(); % run accelerate_simulation.m function
        current_time = current_time + simulation_time + extrapolate_time % update current_time
        XB_iteration_count = XB_iteration_count + 1 % update iteration count
    else % if running another simulation puts current time > original tstop
        run_XB = 0; % tell accelerate_simulation function not to run XBeach
        extrapolate_time = original_tstop - current_time; % set extrapolate_time to difference between current time and original tstop
        accelerate_simulation(); % run accelerate_simulation.m function
        current_time = current_time + extrapolate_time
    end
end

%% Output final profile 

% Import final profile
global extrapolated_profile
z_final = extrapolated_profile; 
if posdwn == 1 % if posdwn is 1
    z_final = -1.*z_final; % flip profile to account for different positive directions
end

figure; % create figure
plot(x, z_initial, 'LineStyle', '--', 'color', 'k', 'linewidth',1.25); % plot initial profile
hold on; % on same figure
plot(x, z_final, 'color', 'b','linewidth',1.25); % plot final profile
xlabel('Cross-shore distance (m)'); ylabel('Elevation (m)'); % add axis labels
legend('Initial Profile', 'Accelerated Simulation', 'Location', 'South'); % add legend

toc; %stop timer