function p_refined = refine_vector(p_input)
    % Refines vector (to double density)
    
    p_refined = zeros(length(p_input)+length(p_input)-1, 1);

    for n = 1:(length(p_input)-1)
        p_refined(2*n-1)  = p_input(n)                     ;
        p_refined(2*n  )  = 0.5*(p_input(n) + p_input(n+1));
    end
    p_refined(end)        = p_input(end)                   ;

    % Test
    % plot(p_input    , p_input   , '+-');
    % hold on;
    % plot(p_refined  , p_refined , 'x-');
end 