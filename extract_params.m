%extract_params
%%% FUNCTION TO EXTRACT DATA FROM PARAMS.TXT FILE %%%#
% Author: Ellie Newell, University of Southampton

function [params_value] = extract_params(variable)
    params = fopen('params.txt'); % open params.txt file
    read_params = textscan(params, '%s %s ', 'Delimiter', '=', 'CommentStyle', '%', 'TextType', 'string'); % read params.txt ignoring % lines 
    fclose(params); % close params.txt
    variable_pat = contains(read_params{1,1}, string(variable)); % array of 1 or 0 if params_data has variable
    in_params = find(variable_pat); % index value of variable location (variable_pat = 1)
    default = ismember(1, variable_pat); % find if variable is in params_data
    if default == 1 % if yes
        params_value = str2double(read_params{1,2}(in_params));  % associated value converted to a double
    elseif default == 0 % if not 
        if strcmp(variable, 'D50') == 1 % if variable is D50
            params_value = 0.0002; % use XBeach default value
        elseif strcmp(variable, 'zs0') == 1 % if variable is zs0
            params_value = 0; % use XBeach default value
        elseif strcmp(variable, 'posdwn') == 1 % if variable is posdwn
            params_value = 1; % use XBeach default value
        elseif strcmp(variable, 'tintg') == 1 % if variable is tintg
            params_value = 1;
        end
    end
end
