% Make coarser mesh
function p_coarsened = coarsen_even_vector(p_input)
    % Coarsens vector (to half density)
    % For vector with even length
    
    % Extract inside points
    p_inside              = p_input(2:end-1);
    % Mid-point them
    p_midpoints           = 0.5 * (p_inside(2:end) + p_inside(1:(end-1)));
    % Coarsen midpoints vector
    p_coarsened_midpoints = coarsen_odd_vector(p_midpoints);
    % Rebuild vector
    p_coarsened     = [p_input(1) ; p_coarsened_midpoints(:) ; p_input(end)];
    % Return with same orientation
    if size(p_input, 1) < size(p_input, 2)
        p_coarsened = transpose(p_coarsened);
    end

    % Test
    % plot(p_input    , p_input    , '+-');
    % hold on;
    % plot(p_coarsened, p_coarsened, 'x-');
end 