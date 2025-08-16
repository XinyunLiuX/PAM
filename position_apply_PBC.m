function wrapped_x = position_apply_PBC(x, box_length)
% Inputs:  
% x = coordinate 
% box_length = size of periodic domain [0, box_length]

% Output: 
% wrapped_x = wrapped x in [0, box_length]

    wrapped_x = mod(x, box_length);
end