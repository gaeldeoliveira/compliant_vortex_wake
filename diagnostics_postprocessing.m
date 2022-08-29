plot(y_range1, u_a1           , '.-'); hold on
plot(y_range1, u_VS_actuator1 , '.-');
plot(y_range1, u_ipc_actuator1, '.-');
grid on;
legend('Free-stream + Wake + Body', 'Wake', 'Body');


% % % For ipc (axi) post-processing
ipc = DSD_axi.ipc
% Get local versors
[n_versor_x, n_versor_y] = normal_versors_column(ipc.aps);

% Define force components on each panel
f_x = ipc.cp .* n_versor_x;
f_y = ipc.cp .* n_versor_y;

% Now integrate using explicit trapeze rule
% Get panel lengths
l_column = panel_lengths_column(ipc.aps);

F_x = sum(f_x .* l_column);
F_y = sum(f_y .* l_column);

% % New stuff (trapezoid integration, prepare for axi)
s_column = cumsum(l_column);
F_x_bis = trapz(s_column, f_x)
F_y_bis = trapz(s_column, f_y)
% % Now try with radial scaling
% Get panel locations (approx)
[x_c_column, y_c_column] = ipc.aps.centerpoint_locations_columns();
% Compute effect of loading 
F_x_bis = trapz(s_column, f_x .* (2 * pi* y_c_column)) ./ trapz(s_column, 2 * pi* 0.6 * ones(size(s_column)))
F_y_bis = trapz(s_column, f_y .* (2 * pi* y_c_column)) ./ trapz(s_column, 2 * pi* 0.6 * ones(size(s_column)))


