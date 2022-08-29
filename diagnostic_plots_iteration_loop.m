% Make verification plots
% u_VS_on_VS   % x velocity of VS on VS points
% v_VS_on_VS   % y velocity of VS on VS points
% u_ipc_on_VS  %
% v_ipc_on_VS  %

% Extract sides
id_top = 1:DSD.n_panels;
id_bot = DSD.n_panels+(1:DSD.n_panels);

% Plot induced velocities ([0 0.4470 0.7410] and [0.8500 0.3250 0.0980])
figure(1)
subplot(211)
plot(DSD.VS.x_center(id_top), u_VS_on_VS(id_top) , 'Color', [0      0.4470 0.7410]); hold on
plot(DSD.VS.x_center(id_bot), u_VS_on_VS(id_bot) , 'Color', [0      0.4470 0.7410]); hold on
plot(DSD.VS.x_center(id_top), u_ipc_on_VS(id_top), 'Color', [0.8500 0.3250 0.0980]); hold on
plot(DSD.VS.x_center(id_bot), u_ipc_on_VS(id_bot), 'Color', [0.8500 0.3250 0.0980]); hold off
subplot(212)
plot(DSD.VS.x_center(id_top), v_VS_on_VS(id_top) , 'Color', [0      0.4470 0.7410]); hold on
plot(DSD.VS.x_center(id_bot), v_VS_on_VS(id_bot) , 'Color', [0      0.4470 0.7410]); hold on
plot(DSD.VS.x_center(id_top), v_ipc_on_VS(id_top), 'Color', [0.8500 0.3250 0.0980]); hold on
plot(DSD.VS.x_center(id_bot), v_ipc_on_VS(id_bot), 'Color', [0.8500 0.3250 0.0980]); hold off
legend('VS_on_VS_top', 'VS_on_VS_bot', 'ipc_on_VS_top', 'ipc_on_VS_bot')

figure(9)


