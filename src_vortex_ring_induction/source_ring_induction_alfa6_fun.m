function [u_z, u_r, u_mag] = source_ring_induction_alfa6_fun(Gamma, Zr, Rr, delta, ellipke_tol, Zt, Rt)
    % Restate induction Point (P) in reference paper notation
    % (embedded in XZ plane, coordinates given as P = [r, 0, z])
    % Gamma       % [m3/s/m ] - (no vectorization)     % Source ring: strenght of source (per unit lenght, may change)
    % Zr          % [m      ] - (no vectorization)     % Source ring: streamwise stance of ring
    % Rr          % [m      ] - (no vectorization)     % Source ring: radius of ring
    % delta       % [adim   ] - (no vectorization)     % Desingularization parameter                (typical: 0  )
    % ellipke_tol % [adim   ] - (no vectorization)     % Toleration of elliptic integral evaluation (typical: eps)
    % Zt          % [m      ] - (vectorizable, tier 1) % Target point: streamwise stance
    % Rt          % [m      ] - (vectorizable, tier 1) % Target point: distane to ring symmetry axis (same as z axis)
    

    % Make radius of target point positive (and store sign to unfold velocity field later)
    % (radius is distance of target point to axis of axi-symmetry of ring)
    Rt0    = Rt;
    Rt     = abs(Rt);
    
    % Shift origin of reference frame to center of ring (axisymmetric
    % assumed, only need to shift along z axis)
    Zt    = Zt - Zr;        % (vectorizable, tier 1)
    
    % Scale target point coordinates to radius of ring
    x_bar = Rt / Rr;        % (vectorizable, tier 1)
    z_bar = Zt / Rr;        % (vectorizable, tier 1) (also denoted as h_bar in deduction)

    % Only proceeed if radius of ring is greater than 0

    if Rr > 0
        % Compute Modulus of Elliptic Integral
        k2       = 4 * x_bar ./ ((1+x_bar).^2 + z_bar.^2);             % (vectorizable, tier 1)
        % Compute Complementary Modulus of Elliptic Integral
        k2_prime = 1 - k2;                                      % (vectorizable, tier 1)
        
        % Compute Elliptic Integrals of the first (K of m) and second (E of m) kind
        % (do this once and reuse result throughout expression, costly operation)
        try
            % Compute two first elliptic integrals:
            [K,E] = ellipke( k2 , ellipke_tol);
        catch
            disp('Error: m values innapropriate');
%             try
%                 disp('Subtract eps(100) to k2 to bypass propagation of float errors')
%                 [K,E] = ellipke( k2 - eps(100) , ellipke_tol);
%             catch
%                 disp('No way forward')
%             end
        end
        % Make third elliptic integral from the two first ones
        D = (K-E) ./ k2;
        
        % Proceed to source ring induction, as deduced by modifying S.J.
        % Newman's deduction for the vortex ring (Southampton TR AFM-11/03)
        
        % Make common multipliers once
        a = Gamma ./ ( pi() .* Rr .* k2_prime);
        b = ((1+x_bar).^2 + z_bar.^2).^(3/2);
        
        % Radial induction % Deprecated: (symmetric of deduction, for consistency with 2d, probably related to change of integration variables (tet to phi, dtet=-dphi) ) (Reason: -gamma effect handled in induction_functions_axi wrapper) 
        u_r = a ./ b    .*   (2 * (K-D) - (1+x_bar) .* E);
        % Axial induction  % Deprecated: (symmetric of deduction, for consistency with 2d, probably related to change of integration variables (tet to phi, dtet=-dphi) ) (Reason: -gamma effect handled in induction_functions_axi wrapper) 
        u_z = a ./ b    .*   (          -    z_bar  .* E);
        % Velocity magnitude
        u_mag = sqrt(u_r.^2 + u_z.^2);
        
        % Restate velocity into (x_hat,y_hat) plane coordinates of external flow analysis
        u_r   = u_r .* sign(Rt0);
    else
        % Ignore
        u_z     = zeros(size(Zt + Rt));
        u_r     = zeros(size(Zt + Rt));
        u_mag   = zeros(size(Zt + Rt));
    end
    
end

%     % Older version, based on a cheap (and incorrect) trick!
%
%     % Rt     % [m       ] - (vectorizable, tier 1) Corresponds to x coordinate of induction point (vectorizable)
%     % Zt     % [m       ] - (vectorizable, tier 1) Corresponds to z coordinate of induction point (vectorizable)
% 
%     % Run vortex ring induction function
%     [u_z_vortex, u_r_vortex, u_mag_vortex] = vortex_ring_induction_alfa6_fun(Gamma, Zr, Rr, delta, ellipke_tol, Zt, Rt);
%     % Switch velocity components
%     u_z =   u_r_vortex .* sign(Rt);
%     u_r = - u_z_vortex .* sign(Rt);
%     % Leave velocity magnitude unaffected
%     u_mag = u_mag_vortex;
%     
%     % The above was wrong!!! Normal direction trick to switch from vortex
%     % to source only works for point sources or specific shapes in planar
%     % flow. 
