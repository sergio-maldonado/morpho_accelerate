%extract_jons
%%% FUNCTION TO EXTRACT WAVE HEIGHT FROM JONSWAP.TXT FILE %%%
% Author: Ellie Newell, University of Southampton

function [params_value] = extract_params(variable)
    jons = fopen('jonswap.txt'); % open jonswap.txt file
    read_jons = textscan(jons, '%s %s ', 'Delimiter', '=', 'CommentStyle', '%', 'TextType', 'string'); % read params.txt ignoring % lines 
    fclose(jons); % jonswap.txt
    variable_pat = contains(read_jons{1,1}, string(variable)); % array of 1 or 0 if jons_data has variable
    in_params = find(variable_pat); % index value of variable location (variable_pat = 1)
    default = ismember(1, variable_pat); % find if variable is in jons_data
    params_value = str2double(read_jons{1,2}(in_params));  % associated value converted to a double
end
