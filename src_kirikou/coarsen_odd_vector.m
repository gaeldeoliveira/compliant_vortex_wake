function p_coarsened = coarsen_odd_vector(p_input)
    % Coarsens vector (to half density)
    % For vector with odd length
    p_coarsened = p_input(1:2:length(p_input));

end