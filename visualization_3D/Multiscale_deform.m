% ---------------------------------------------------------------------------------- %
% Multiscale structure deformation visualization with microstructure overlay.
% ---------------------------------------------------------------------------------- %
%% Define structure and loading parameters
nx=40; ny=30; nz=40; nelx=20; nely=20; nelz=20;
Tt=2;  Tb=3;
x=zeros(ny,nx,nz); x([1:Tt,end-(Tb-1):end],:,:)=1;  
load('3Dsupported.mat','Phi0'); Phi=reshape(Phi0,[nely+1,nelx+1,nelz+1]);
load('U.mat','U');
%% Color map
color=[0 0 223;     0 0 239;     0 0 255;     0 16 255;    0 32 255;   0 48 255;   0 64 255;   0 80 255;    0 96 255;  
        0 111 255;   0 128 255;   0 143 255;   0 159 255;   0 175 255;  0 191 255;  0 207 255;  0 223 255;   0 239 255;   
        0 255 255;   16 255 239;  32 255 223;  48 255 207;  64 255 191; 80 255 175; 96 255 159; 111 255 143; 128 255 128; 
        143 255 111; 159 255 96;  175 255 80;  191 255 64;  207 255 48; 223 255 32; 239 255 16; 255 239 0;   255 223 0;   
        255 207 0;   255 191 0;   255 175 0;   255 159 0;   255 143 0;  255 128 0;  255 111 0;  255 96 0;    255 80 0;   
        255 64 0;    255 48 0;    255 32 0;    255 16 0;    255 0 0;    239 0 0;    223 0 0;    207 0 0;     191 0 0;]./255;
figure(1); clf; deformation_multiscale(x,Phi,U,color);

%% Core multiscale deformation visualization
function eigenmode_disp_3D(rho,U,A0,color)
[nely,nelx,nelz] = size(rho);
hx = 1; hy = 1; hz = 1;            % User-defined unit element size
face = [1 2 3 4; 2 6 7 3; 4 3 7 8; 1 5 8 4; 1 2 6 5; 5 6 7 8];
% set(gcf,'Name','ISO display','NumberTitle','off');
for k = 1:nelz
    z = (k-1)*hz;
    for i = 1:nelx
        x = (i-1)*hx;
        for j = 1:nely
            y = nely*hy - (j-1)*hy;
            n1 = (nelx+1)*(nely+1)*(k-1) + (nely+1)*(i-1) + (j+1);  %(02)
            n2 = (nelx+1)*(nely+1)*(k-1) + (nely+1)*(i)   + (j+1);  %(04)
            n3 = (nelx+1)*(nely+1)*(k-1) + (nely+1)*(i)   + j;      %(03)
            n4 = (nelx+1)*(nely+1)*(k-1) + (nely+1)*(i-1) + j;      %(01)
            n5 = (nelx+1)*(nely+1)*(k)   + (nely+1)*(i-1) + (j+1);  %(12)
            n6 = (nelx+1)*(nely+1)*(k)   + (nely+1)*(i)   + (j+1);  %(14)
            n7 = (nelx+1)*(nely+1)*(k)   + (nely+1)*(i)   + j;      %(13)
            n8 = (nelx+1)*(nely+1)*(k)   + (nely+1)*(i-1) + j;      %(11)
            % edof = [3*n1-2 3*n1-1 3*n1; 3*n2-2 3*n2-1 3*n2;
            %         3*n3-2 3*n3-1 3*n3; 3*n4-2 3*n4-1 3*n4;
            %         3*n5-2 3*n5-1 3*n5; 3*n6-2 3*n6-1 3*n6;
            %         3*n7-2 3*n7-1 3*n7; 3*n8-2 3*n8-1 3*n8];
            edof = [3*n4-2 3*n4-1 3*n4; 3*n1-2 3*n1-1 3*n1; 
                    3*n2-2 3*n2-1 3*n2; 3*n3-2 3*n3-1 3*n3; 
                    3*n8-2 3*n8-1 3*n8; 3*n5-2 3*n5-1 3*n5;
                    3*n6-2 3*n6-1 3*n6; 3*n7-2 3*n7-1 3*n7; ];
            edof = edof(:,[3 1 2]);
            Ue = U(edof); Ue_node=vecnorm(Ue,2,2);
            if (rho(j,i,k) > 0.01)  % User-defined display density threshold
                vert = [x y z; x y-hx z; x+hx y-hx z; x+hx y z; 
                        x y z+hx;x y-hx z+hx; x+hx y-hx z+hx;x+hx y z+hx];
                vert(:,[2 3]) = vert(:,[3 2]); %vert(:,2,:) = -vert(:,2,:);
                vert = vert+Ue/2;
                if rho(j,i,k)==1; RGB=[191/255 0 0];transp=1; else; RGB=color(ceil(A0*54),:);transp=0.6; end
                patch('Faces',face,'Vertices',vert,'FaceColor',RGB,'edgecolor','none');
                if x==0||x==nelx-1||y==nely||y==1||z==0||z==nelz-1
                patch('Faces',face,'Vertices',vert,'FaceVertexCData',Ue_node,'FaceColor','interp','edgecolor','none','facealpha',0.55);
                end
                hold on;
            end
        end
    end
end
axis equal; axis tight; axis off; box on; view([30,30]); pause(1e-6);
end

% function deformation_multiscale(den,Phi,U,color)
% [nely,nelx] = size(den);
% hx = 1; hy = 1;
% face = [1 2 3 4];
% scale_def = 30;
% 
% % Downsample microstructure for visualization
% stride = 1; % downsample per stride
% Phi_v = Phi(1:stride:end, 1:stride:end);
% [nely_phi, nelx_phi] = size(Phi_v);
% nely_m = nely_phi - 1;
% nelx_m = nelx_phi - 1;
% 
% % Precompute micro-cell solid flags and local coords
% phi_cells = 0.25 * (Phi_v(1:end-1,1:end-1) + Phi_v(2:end,1:end-1) + ...
%                     Phi_v(1:end-1,2:end) + Phi_v(2:end,2:end));
% [solid_rows, solid_cols] = find(phi_cells > 0);
% num_solid = length(solid_rows);
% 
% if num_solid == 0
%     return;
% end
% 
% % Precompute local (xi,eta) corners for each micro cell (0 to 1 range)
% xi_l = (solid_cols-1) / nelx_m;
% xi_r = solid_cols / nelx_m;
% eta_b = 1 - solid_rows / nely_m;
% eta_t = 1 - (solid_rows-1) / nely_m;
% 
% hold on;
% for i = 1:nelx
%     for j = 1:nely
% 
%         % Macro element node indices
%         n1 = (nely+1)*(i-1) + j;
%         n2 = (nely+1)*i + j;
%         Ue = [U(2*n1-1) U(2*n1); U(2*n1+1) U(2*n1+2); ...
%               U(2*n2+1) U(2*n2+2); U(2*n2-1) U(2*n2)];
% 
%         % Macro element undeformed corners: [BL; TL; TR; BR]
%         xL = (i-1)*hx; xR = i*hx;
%         yB = nely*hy - j*hy; yT = nely*hy - (j-1)*hy;
%         P = [xL yT; xL yB; xR yB; xR yT] + Ue * scale_def;
% 
%         % --- Microstructure patches ---
%         all_verts = zeros(num_solid*4, 2);
%         all_cdata = zeros(num_solid, 1);
% 
%         % Local cell corner xi,eta arrays
%         xi_quad = [xi_l, xi_r, xi_r, xi_l];
%         eta_quad = [eta_t, eta_t, eta_b, eta_b];
% 
%         for k = 1:num_solid
%             xi = xi_quad(k,:);
%             eta = eta_quad(k,:);
%             % Bilinear interpolation of 4 corners
%             N1 = (1-xi).*(1-eta); N2 = (1-xi).*eta;
%             N3 = xi.*eta;         N4 = xi.*(1-eta);
%             vert_def = [N1*P(1,1)+N2*P(2,1)+N3*P(3,1)+N4*P(4,1);
%                         N1*P(1,2)+N2*P(2,2)+N3*P(3,2)+N4*P(4,2)]';
%             all_verts((k-1)*4+1:k*4, :) = vert_def;
% 
%             % Color from displacement magnitude at cell center
%             xi_c = (xi_l(k)+xi_r(k))/2;
%             eta_c = (eta_t(k)+eta_b(k))/2;
%             N1c = (1-xi_c)*(1-eta_c); N2c = (1-xi_c)*eta_c;
%             N3c = xi_c*eta_c;         N4c = xi_c*(1-eta_c);
%             u_cen = N1c*Ue(1,:) + N2c*Ue(2,:) + N3c*Ue(3,:) + N4c*Ue(4,:);
%             all_cdata(k) = norm(u_cen);
%         end
%         if den(j,i)==0
%         faces = reshape(1:4*num_solid, 4, num_solid)';
%         patch('Vertices', all_verts, 'Faces', faces, ...
%               'FaceVertexCData', all_cdata, ...
%               'FaceColor', 'flat', 'EdgeColor', 'none');
%         else
%         % --- Macro element translucent overlay ---
%         Ue_node = vecnorm(Ue,2,2);
%         patch('Faces', face, 'Vertices', P, ...
%               'FaceVertexCData', Ue_node, 'FaceColor', 'interp', ...
%               'EdgeColor', 'none');
%         end
%     end
% end
% 
% axis equal; axis tight; axis off; box off; grid off;
% colormap(color);
% end