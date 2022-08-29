
% % % Test induction functions
N_range           = 200;
x_range           = linspace(-4 , 4, N_range);
y_range           = linspace(-4 , 4, N_range);
[x_mesh , y_mesh]= meshgrid(x_range, y_range);

delta             = 0;
ellipke_tol       = eps();

Gamma             = 1;
Zr                = 0;
Rr                = 0.5;

u_z_mesh          = zeros(size(x_range)); 
u_r_mesh          = zeros(size(y_range));
u_mag_mesh        = zeros(size(y_range));

% Sequential Evaluation over mesh of targets
for i = 1:N_range
    for j = 1:N_range
        [u_z_mesh(i,j), u_r_mesh(i,j), u_mag_mesh(i,j)] = source_ring_induction_alfa6_fun(Gamma, Zr, Rr, delta, ellipke_tol, x_mesh(i,j), y_mesh(i,j));
    end
end
% Vectorized Evaluation over mesh of targets (now via wrapper
[u_z_meshB, u_r_meshB, u_mag_meshB] = source_ring_induction_alfa6_fun(Gamma, Zr, Rr, delta, ellipke_tol, x_mesh, y_mesh);
% This should be empty, when all works well with the vectorization
find(not(u_z_mesh == -u_z_meshB))
find(not(u_r_mesh == -u_r_meshB))
find(not(u_mag_mesh == u_mag_meshB))
% Wrapper test
[u_z_meshC, u_r_meshC] = induction_functions_axi.constant_strenght_singular_source_ring_axi(x_mesh, y_mesh, Zr, Rr, Gamma);
% This should be empty, when all works well with the wrapper
find(not(u_z_mesh == -u_z_meshC))
find(not(u_r_mesh == -u_r_meshC))



figure(1)
surf(x_mesh, y_mesh, u_z_meshB); view(2); shading flat
figure(2)
surf(x_mesh, y_mesh, u_r_meshB); view(2); shading flat
figure(3)
quiver(x_mesh, y_mesh, u_z_meshB, u_r_meshB);
axis([-0.2 0.2 -1 1])
figure(4)
tet_range = linspace(0,2*pi(), 36); streamline(x_mesh, y_mesh, u_z_meshB, u_r_meshB, 4*cos(tet_range), 4*sin(tet_range))
figure(5)
contour(x_mesh, y_mesh, u_mag_meshB)

[u_z_meshV, u_r_meshV     ] = induction_functions.constant_strenght_singular_source_pair(x_mesh, y_mesh, Zr, Rr, Gamma);
figure(11)
surf(x_mesh, y_mesh, u_z_meshV); view(2); shading flat
figure(21)
surf(x_mesh, y_mesh, u_r_meshV); view(2); shading flat
figure(31)
quiver(x_mesh, y_mesh, u_z_meshV, u_r_meshV);
axis([-0.2 0.2 -1 1])


% % Correspondence to panel method
x and u_r (axi-symmetric integration) goes to y and u_y or v (panel)
z and u_z (axi-symmetric integration) goes to x and u_x or u (panel)