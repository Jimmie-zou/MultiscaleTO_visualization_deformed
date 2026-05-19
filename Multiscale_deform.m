% ---------------------------------------------------------------------------------- %
% Multiscale structure deformation visualization with microstructure overlay.
% ---------------------------------------------------------------------------------- %
%% Define structure and loading parameters
nx=50; ny=30; nelx=50; nely=50;
Tt=2;  Tb=4;
x=zeros(ny,nx); x([1:Tt,end-(Tb-1):end],:)=1; 
load('lattice.mat','Phi');
load('U.mat','U');
Phi = reshape(Phi,nely+1,nelx+1);
%% Color map
color=[0 0 223;     0 0 239;     0 0 255;     0 16 255;    0 32 255;   0 48 255;   0 64 255;   0 80 255;    0 96 255;  
        0 111 255;   0 128 255;   0 143 255;   0 159 255;   0 175 255;  0 191 255;  0 207 255;  0 223 255;   0 239 255;   
        0 255 255;   16 255 239;  32 255 223;  48 255 207;  64 255 191; 80 255 175; 96 255 159; 111 255 143; 128 255 128; 
        143 255 111; 159 255 96;  175 255 80;  191 255 64;  207 255 48; 223 255 32; 239 255 16; 255 239 0;   255 223 0;   
        255 207 0;   255 191 0;   255 175 0;   255 159 0;   255 143 0;  255 128 0;  255 111 0;  255 96 0;    255 80 0;   
        255 64 0;    255 48 0;    255 32 0;    255 16 0;    255 0 0;    239 0 0;    223 0 0;    207 0 0;     191 0 0;]./255;
figure(1); clf; deformation_multiscale(x,Phi,U,color);


%% Core multiscale deformation visualization
function deformation_multiscale(den,Phi,U,color)
[nely,nelx] = size(den);
hx = 1; hy = 1;
face = [1 2 3 4];
scale_def = 30;

% Downsample microstructure for visualization
stride = 1; % downsample per stride
Phi_v = Phi(1:stride:end, 1:stride:end);
[nely_phi, nelx_phi] = size(Phi_v);
nely_m = nely_phi - 1;
nelx_m = nelx_phi - 1;

% Precompute micro-cell solid flags and local coords
phi_cells = 0.25 * (Phi_v(1:end-1,1:end-1) + Phi_v(2:end,1:end-1) + ...
                    Phi_v(1:end-1,2:end) + Phi_v(2:end,2:end));
[solid_rows, solid_cols] = find(phi_cells > 0);
num_solid = length(solid_rows);

if num_solid == 0
    return;
end

% Precompute local (xi,eta) corners for each micro cell (0 to 1 range)
xi_l = (solid_cols-1) / nelx_m;
xi_r = solid_cols / nelx_m;
eta_b = 1 - solid_rows / nely_m;
eta_t = 1 - (solid_rows-1) / nely_m;

hold on;
for i = 1:nelx
    for j = 1:nely

        % Macro element node indices
        n1 = (nely+1)*(i-1) + j;
        n2 = (nely+1)*i + j;
        Ue = [U(2*n1-1) U(2*n1); U(2*n1+1) U(2*n1+2); ...
              U(2*n2+1) U(2*n2+2); U(2*n2-1) U(2*n2)];

        % Macro element undeformed corners: [BL; TL; TR; BR]
        xL = (i-1)*hx; xR = i*hx;
        yB = nely*hy - j*hy; yT = nely*hy - (j-1)*hy;
        P = [xL yT; xL yB; xR yB; xR yT] + Ue * scale_def;

        % --- Microstructure patches ---
        all_verts = zeros(num_solid*4, 2);
        all_cdata = zeros(num_solid, 1);

        % Local cell corner xi,eta arrays
        xi_quad = [xi_l, xi_r, xi_r, xi_l];
        eta_quad = [eta_t, eta_t, eta_b, eta_b];

        for k = 1:num_solid
            xi = xi_quad(k,:);
            eta = eta_quad(k,:);
            % Bilinear interpolation of 4 corners
            N1 = (1-xi).*(1-eta); N2 = (1-xi).*eta;
            N3 = xi.*eta;         N4 = xi.*(1-eta);
            vert_def = [N1*P(1,1)+N2*P(2,1)+N3*P(3,1)+N4*P(4,1);
                        N1*P(1,2)+N2*P(2,2)+N3*P(3,2)+N4*P(4,2)]';
            all_verts((k-1)*4+1:k*4, :) = vert_def;

            % Color from displacement magnitude at cell center
            xi_c = (xi_l(k)+xi_r(k))/2;
            eta_c = (eta_t(k)+eta_b(k))/2;
            N1c = (1-xi_c)*(1-eta_c); N2c = (1-xi_c)*eta_c;
            N3c = xi_c*eta_c;         N4c = xi_c*(1-eta_c);
            u_cen = N1c*Ue(1,:) + N2c*Ue(2,:) + N3c*Ue(3,:) + N4c*Ue(4,:);
            all_cdata(k) = norm(u_cen);
        end
        if den(j,i)==0
        faces = reshape(1:4*num_solid, 4, num_solid)';
        patch('Vertices', all_verts, 'Faces', faces, ...
              'FaceVertexCData', all_cdata, ...
              'FaceColor', 'flat', 'EdgeColor', 'none');
        else
        % --- Macro element translucent overlay ---
        Ue_node = vecnorm(Ue,2,2);
        patch('Faces', face, 'Vertices', P, ...
              'FaceVertexCData', Ue_node, 'FaceColor', 'interp', ...
              'EdgeColor', 'none');
        end
    end
end

axis equal; axis tight; axis off; box off; grid off;
colormap(color);
end