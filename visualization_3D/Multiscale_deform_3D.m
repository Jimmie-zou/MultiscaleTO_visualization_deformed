% ---------------------------------------------------------------------------------- %
% 3D Multiscale structure deformation visualization with microstructure overlay.
% ---------------------------------------------------------------------------------- %
%% Define structure and loading parameters
nx=40; ny=30; nz=40; nelx=20; nely=20; nelz=20;
Tt=2;  Tb=3;
x=zeros(ny,nx,nz); x([1:Tt,end-(Tb-1):end],:,:)=1;
load('3Dsupported.mat','Phi0'); 
Phi=reshape(Phi0,[nely+1,nelx+1,nelz+1]); Phi=permute(Phi,[3 2 1]);
load('U.mat','U'); U = 10*U;

%% Color map
color=[0 0 223;     0 0 239;     0 0 255;     0 16 255;    0 32 255;   0 48 255;   0 64 255;   0 80 255;    0 96 255;
        0 111 255;   0 128 255;   0 143 255;   0 159 255;   0 175 255;  0 191 255;  0 207 255;  0 223 255;   0 239 255;
        0 255 255;   16 255 239;  32 255 223;  48 255 207;  64 255 191; 80 255 175; 96 255 159; 111 255 143; 128 255 128;
        143 255 111; 159 255 96;  175 255 80;  191 255 64;  207 255 48; 223 255 32; 239 255 16; 255 239 0;   255 223 0;
        255 207 0;   255 191 0;   255 175 0;   255 159 0;   255 143 0;  255 128 0;  255 111 0;  255 96 0;    255 80 0;
        255 64 0;    255 48 0;    255 32 0;    255 16 0;    255 0 0;    239 0 0;    223 0 0;    207 0 0;     191 0 0]./255;
figure(1); clf;
deformation_multiscale_3D(x,Phi,U,color);

%% Core 3D multiscale deformation visualization (optimized)
function deformation_multiscale_3D(den,Phi,U,color)
[nely, nelx, nelz] = size(den);
hx = 1; hy = 1; hz = 1;
scale_def = 30;

% --- Microstructure level-set on fine grid ---
stride = 2;
Phi_v = Phi(1:stride:end, 1:stride:end, 1:stride:end);
[nely_phi, nelx_phi, nelz_phi] = size(Phi_v);
nely_m = nely_phi - 1; nelx_m = nelx_phi - 1; nelz_m = nelz_phi - 1;

phi_cells = (Phi_v(1:end-1,1:end-1,1:end-1) + Phi_v(2:end,1:end-1,1:end-1) + ...
             Phi_v(1:end-1,2:end,1:end-1) + Phi_v(2:end,2:end,1:end-1) + ...
             Phi_v(1:end-1,1:end-1,2:end) + Phi_v(2:end,1:end-1,2:end) + ...
             Phi_v(1:end-1,2:end,2:end) + Phi_v(2:end,2:end,2:end)) / 8;

[solid_y, solid_x, solid_z] = ind2sub(size(phi_cells), find(phi_cells > 0));
num_solid = length(solid_y);
if num_solid == 0
    return;
end

% --- Micro cell local coordinates within macro element [0,1]^3 ---
xi_l   = (solid_x-1) / nelx_m;    xi_r   = solid_x     / nelx_m;
eta_b  = (solid_y-1) / nely_m;    eta_t  = solid_y     / nely_m;
zeta_b = (solid_z-1) / nelz_m;    zeta_f = solid_z     / nelz_m;

% --- Pre-compute shape functions (independent of macro element) ---
xi_8   = [xi_l, xi_r, xi_r, xi_l, xi_l, xi_r, xi_r, xi_l];
eta_8  = [eta_b, eta_b, eta_t, eta_t, eta_b, eta_b, eta_t, eta_t];
zeta_8 = [zeta_b, zeta_b, zeta_b, zeta_b, zeta_f, zeta_f, zeta_f, zeta_f];

N1 = (1-xi_8).*(1-eta_8).*(1-zeta_8);
N2 = xi_8.*(1-eta_8).*(1-zeta_8);
N3 = xi_8.*eta_8.*(1-zeta_8);
N4 = (1-xi_8).*eta_8.*(1-zeta_8);
N5 = (1-xi_8).*(1-eta_8).*zeta_8;
N6 = xi_8.*(1-eta_8).*zeta_8;
N7 = xi_8.*eta_8.*zeta_8;
N8 = (1-xi_8).*eta_8.*zeta_8;

xi_cen = (xi_l + xi_r) / 2;
eta_cen = (eta_b + eta_t) / 2;
zeta_cen = (zeta_b + zeta_f) / 2;

N1c = (1-xi_cen).*(1-eta_cen).*(1-zeta_cen);
N2c = xi_cen.*(1-eta_cen).*(1-zeta_cen);
N3c = xi_cen.*eta_cen.*(1-zeta_cen);
N4c = (1-xi_cen).*eta_cen.*(1-zeta_cen);
N5c = (1-xi_cen).*(1-eta_cen).*zeta_cen;
N6c = xi_cen.*(1-eta_cen).*zeta_cen;
N7c = xi_cen.*eta_cen.*zeta_cen;
N8c = (1-xi_cen).*eta_cen.*zeta_cen;

% --- Pre-compute face connectivity (vectorized, no loop) ---
face_local = [1 2 3 4; 5 6 7 8; 1 2 6 5; 2 3 7 6; 3 4 8 7; 4 1 5 8];
micro_offsets = (0:num_solid-1)' * 8;
solid_faces = repmat(face_local, num_solid, 1) + ...
              reshape(repmat(micro_offsets, 1, 6)', [], 1);

% --- Batching: reduce patch() calls ---
VOID_BATCH = 200;
SOLID_BATCH = 1000;

n_void = 0; v_off = 0;
void_v = cell(VOID_BATCH, 1);
void_f = cell(VOID_BATCH, 1);
void_c = cell(VOID_BATCH, 1);

n_solid = 0; s_off = 0;
solid_v = cell(SOLID_BATCH, 1);
solid_f = cell(SOLID_BATCH, 1);
solid_c = cell(SOLID_BATCH, 1);

flush_void = @(vv, ff, cc, n) patch('Vertices', cat(1, vv{1:n}), ...
    'Faces', cat(1, ff{1:n}), 'FaceVertexCData', cat(1, cc{1:n}), ...
    'FaceColor', 'flat', 'EdgeColor', 'none');
flush_solid = @(vv, ff, cc, n) patch('Faces', cat(1, ff{1:n}), ...
    'Vertices', cat(1, vv{1:n}), 'FaceVertexCData', cat(1, cc{1:n}), ...
    'FaceColor', 'interp', 'EdgeColor', 'none');

hold on;
for i = 1:nelx
    for j = 1:nely
        for k = 1:nelz
            if i == 1|| j==1||k==1
            % --- Macro element node indices ---
            n1 = (k-1)*(nelx+1)*(nely+1) + (i-1)*(nely+1) + j;
            n2 = (k-1)*(nelx+1)*(nely+1) + i*(nely+1) + j;
            n3 = (k-1)*(nelx+1)*(nely+1) + i*(nely+1) + j+1;
            n4 = (k-1)*(nelx+1)*(nely+1) + (i-1)*(nely+1) + j+1;
            n5 = k*(nelx+1)*(nely+1) + (i-1)*(nely+1) + j;
            n6 = k*(nelx+1)*(nely+1) + i*(nely+1) + j;
            n7 = k*(nelx+1)*(nely+1) + i*(nely+1) + j+1;
            n8 = k*(nelx+1)*(nely+1) + (i-1)*(nely+1) + j+1;

            Ue = [U(3*n1-2) U(3*n1-1) U(3*n1);
                  U(3*n2-2) U(3*n2-1) U(3*n2);
                  U(3*n3-2) U(3*n3-1) U(3*n3);
                  U(3*n4-2) U(3*n4-1) U(3*n4);
                  U(3*n5-2) U(3*n5-1) U(3*n5);
                  U(3*n6-2) U(3*n6-1) U(3*n6);
                  U(3*n7-2) U(3*n7-1) U(3*n7);
                  U(3*n8-2) U(3*n8-1) U(3*n8)];

            P = [(i-1)*hx (j-1)*hy (k-1)*hz; i*hx (j-1)*hy (k-1)*hz;
                 i*hx j*hy (k-1)*hz; (i-1)*hx j*hy (k-1)*hz;
                 (i-1)*hx (j-1)*hy k*hz; i*hx (j-1)*hy k*hz;
                 i*hx j*hy k*hz; (i-1)*hx j*hy k*hz] + Ue * scale_def;

            if den(j,i,k) == 0
                % --- Void element: draw microstructure ---
                verts_x = N1*P(1,1) + N2*P(2,1) + N3*P(3,1) + N4*P(4,1) + ...
                          N5*P(5,1) + N6*P(6,1) + N7*P(7,1) + N8*P(8,1);
                verts_y = N1*P(1,2) + N2*P(2,2) + N3*P(3,2) + N4*P(4,2) + ...
                          N5*P(5,2) + N6*P(6,2) + N7*P(7,2) + N8*P(8,2);
                verts_z = N1*P(1,3) + N2*P(2,3) + N3*P(3,3) + N4*P(4,3) + ...
                          N5*P(5,3) + N6*P(6,3) + N7*P(7,3) + N8*P(8,3);

                u_cen_x = N1c*Ue(1,1) + N2c*Ue(2,1) + N3c*Ue(3,1) + N4c*Ue(4,1) + ...
                          N5c*Ue(5,1) + N6c*Ue(6,1) + N7c*Ue(7,1) + N8c*Ue(8,1);
                u_cen_y = N1c*Ue(1,2) + N2c*Ue(2,2) + N3c*Ue(3,2) + N4c*Ue(4,2) + ...
                          N5c*Ue(5,2) + N6c*Ue(6,2) + N7c*Ue(7,2) + N8c*Ue(8,2);
                u_cen_z = N1c*Ue(1,3) + N2c*Ue(2,3) + N3c*Ue(3,3) + N4c*Ue(4,3) + ...
                          N5c*Ue(5,3) + N6c*Ue(6,3) + N7c*Ue(7,3) + N8c*Ue(8,3);

                n_void = n_void + 1;
                void_v{n_void} = [reshape(verts_x', [], 1), reshape(verts_y', [], 1), reshape(verts_z', [], 1)];
                void_f{n_void} = solid_faces + v_off;
                void_c{n_void} = repelem(sqrt(u_cen_x.^2 + u_cen_y.^2 + u_cen_z.^2), 8);
                v_off = v_off + size(void_v{n_void}, 1);

                if n_void == VOID_BATCH
                    flush_void(void_v, void_f, void_c, VOID_BATCH);
                    n_void = 0; v_off = 0;
                end

            else
                % --- Solid element: skip internal faces ---
                vis = true(6, 1);
                if i > 1  && den(j,i-1,k) ~= 0, vis(6) = false; end
                if i < nelx && den(j,i+1,k) ~= 0, vis(4) = false; end
                if j > 1  && den(j-1,i,k) ~= 0, vis(3) = false; end
                if j < nely && den(j+1,i,k) ~= 0, vis(5) = false; end
                if k > 1  && den(j,i,k-1) ~= 0, vis(1) = false; end
                if k < nelz && den(j,i,k+1) ~= 0, vis(2) = false; end

                vfaces = face_local(vis, :);
                if ~isempty(vfaces)
                    n_solid = n_solid + 1;
                    solid_v{n_solid} = P;
                    solid_f{n_solid} = vfaces + s_off;
                    solid_c{n_solid} = vecnorm(Ue, 2, 2);
                    s_off = s_off + 8;

                    if n_solid == SOLID_BATCH
                        flush_solid(solid_v, solid_f, solid_c, SOLID_BATCH);
                        n_solid = 0; s_off = 0;
                    end
                end
            end
            end
        end
    end
end

% Flush remaining elements
if n_void > 0
    flush_void(void_v, void_f, void_c, n_void);
end
if n_solid > 0
    flush_solid(solid_v, solid_f, solid_c, n_solid);
end

axis equal; axis tight; axis off; box off; grid off; view([-10,-60]); pause(1e-6);
colormap(color);
camlight; lighting flat;
end
