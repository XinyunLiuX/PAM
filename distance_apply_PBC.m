function wrapped_delta = distance_apply_PBC(delta, box_length)
% Inputs: 
% delta = x' - x 
% box_length = size of periodic domain [0, box_length]


% Output: 
% wrapped_delta = wrapped_delta, choosing the smallest distance in periodic 
% boundary condition, including both negative and positive values

    wrapped_delta = delta - box_length*round(delta/box_length);
end