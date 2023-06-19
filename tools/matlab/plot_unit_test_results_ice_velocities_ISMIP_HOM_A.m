clc
clear all
close all

foldername = '../../results_ISMIP_HOM_A';

%% Velocity maps

Ls    = [160,80,40,20,10,5];
umins = [  0,  0, 0, 0, 0, 0];
umaxs = [125,100,80,60,40,30];

for li = 1: length( Ls)
  
  L    = Ls(    li);
  umin = umins( li);
  umax = umaxs( li);

  filename = [foldername '/test_ISMIP_HOM_A_' num2str( L) 'km_output.nc'];

  mesh = read_mesh_from_file( filename);

  %% Plot

  wa = 300;
  ha = 300;
  margins_hor = [25, 25, 25, 100];
  margins_ver = [50, 25];
  H = setup_multipanel_figure( wa, ha, margins_hor, margins_ver);

  for i = 1: size( H.Ax,1)
    for j = 1: size( H.Ax,2)
      set( H.Ax{ 1,j},'xtick',[],'ytick',[],'clim',[umin,umax],'fontsize',24);
    end
  end

  title( H.Ax{ 1,1},'SIA/SSA');
  title( H.Ax{ 1,2},'DIVA');
  title( H.Ax{ 1,3},'BPA');
  pos = get( H.Ax{ 1,3},'position');
  H.Cbar = colorbar( H.Ax{ 1,3},'location','eastoutside');
  ylabel( H.Cbar,'Surface velocity (m/yr)')
  set( H.Ax{ 1,3},'position',pos);

  % SIA/SSA
  u_3D_b_SIASSA = ncread( filename,'u_3D_b_SIASSA');
  plot_mesh_data_b( H.Ax{ 1,1}, mesh, u_3D_b_SIASSA( :,1));

  % DIVA
  u_3D_b_DIVA = ncread( filename,'u_3D_b_DIVA');
  plot_mesh_data_b( H.Ax{ 1,2}, mesh, u_3D_b_DIVA( :,1));

  % BPA
  u_3D_b_BPA = ncread( filename,'u_3D_b_BPA');
  plot_mesh_data_b( H.Ax{ 1,3}, mesh, u_3D_b_BPA( :,1));

end

%% Transects, comparing to the Pattyn et al. (2008) ensemble

% Read Pattyn et al. (2008) data
foldername_ensemble ='/Users/berends/Documents/Datasets/ISMIP-HOM/tc-2-95-2008-supplement/ismip_all/';
[R,model_types] = read_Pattyn2008_ensemble( foldername_ensemble);

%% Plot Pattyn et al. (2008) data
wa = 400;
ha = 300;
margins_hor = [90,75,75,25];
margins_ver = [50,60,70];
H2 = setup_multipanel_figure( wa, ha, margins_hor, margins_ver);

for i = 1: size( H2.Ax,1)
  for j = 1: size( H2.Ax,2)
    set( H2.Ax{ i,j},'fontsize',24,'xlim',[0,1],'xgrid','on','ygrid','on');
  end
end

title( H2.Ax{ 1,1},'160 km')
title( H2.Ax{ 1,2},'80 km')
title( H2.Ax{ 1,3},'40 km')
title( H2.Ax{ 2,1},'20 km')
title( H2.Ax{ 2,2},'10 km')
title( H2.Ax{ 2,3},'5 km')

ylabel( H2.Ax{ 1,1},'Surface x-velocity (m/yr)')
ylabel( H2.Ax{ 2,1},'Surface x-velocity (m/yr)')

xlabel( H2.Ax{ 2,1},'x / L')
xlabel( H2.Ax{ 2,2},'x / L')
xlabel( H2.Ax{ 2,3},'x / L')

set( H2.Ax{ 1,1},'ylim',[umins(1),umaxs(1)]);
set( H2.Ax{ 1,2},'ylim',[umins(2),umaxs(2)]);
set( H2.Ax{ 1,3},'ylim',[umins(3),umaxs(3)]);
set( H2.Ax{ 2,1},'ylim',[umins(4),umaxs(4)]);
set( H2.Ax{ 2,2},'ylim',[umins(5),umaxs(5)]);
set( H2.Ax{ 2,3},'ylim',[umins(6),umaxs(6)]);

% Legend
c_FS     = [0.2,0.5,1.0];
c_HO     = [0.1,0.7,0.3];
c_SIASSA = [1.0,0.0,0.0];
c_DIVA   = [1.0,0.7,0.0];
c_BPA    = [0.5,0.0,0.7];
line(  H2.Ax{1,1},'xdata'   ,[],'ydata',[],'color'    ,c_HO    ,'linewidth',3);
patch( H2.Ax{1,1},'vertices',[],'faces',[],'facecolor',c_HO    ,'edgecolor','none','facealpha',0.7)
line(  H2.Ax{1,1},'xdata'   ,[],'ydata',[],'color'    ,c_FS    ,'linewidth',3);
patch( H2.Ax{1,1},'vertices',[],'faces',[],'facecolor',c_FS    ,'edgecolor','none','facealpha',0.7)
line(  H2.Ax{1,1},'xdata'   ,[],'ydata',[],'color'    ,c_SIASSA,'linewidth',3);
line(  H2.Ax{1,1},'xdata'   ,[],'ydata',[],'color'    ,c_DIVA  ,'linewidth',3);
line(  H2.Ax{1,1},'xdata'   ,[],'ydata',[],'color'    ,c_BPA   ,'linewidth',3);

% Results
for xi = 1:6
  
  if xi == 1
    ex = 'a160';
    i = 1;
    j = 1;
  elseif xi == 2
    ex = 'a080';
    i = 1;
    j = 2;
  elseif xi == 3
    ex = 'a040';
    i = 1;
    j = 3;
  elseif xi == 4
    ex = 'a020';
    i = 2;
    j = 1;
  elseif xi == 5
    ex = 'a010';
    i = 2;
    j = 2;
  elseif xi == 6
    ex = 'a005';
    i = 2;
    j = 3;
  end
  
  x_FS = linspace(0,1,200)';
  u_FS = [];
  u_HO = [];
  
  patch_HO = patch(H2.Ax{i,j},'xdata',[],'ydata',[],'facecolor',c_HO,'edgecolor','none','facealpha',0.5);
  patch_FS = patch(H2.Ax{i,j},'xdata',[],'ydata',[],'facecolor',c_FS,'edgecolor','none','facealpha',0.7);
  line_HO  = line('parent',H2.Ax{i,j},'xdata',[],'ydata',[],'color',c_HO,'linewidth',3);
  line_FS  = line('parent',H2.Ax{i,j},'xdata',[],'ydata',[],'color',c_FS,'linewidth',3);
  
  flds = fields(R.(ex));
  for mi = 1:length(flds)
    m = flds{mi};
    
    % Get transect at y = 0.25
    x = R.(ex).(m).x(:,1);
    y = R.(ex).(m).y(1,:)';
    yi = round(0.25*length(y));
    u = R.(ex).(m).u(:,yi);

    % Determine if this model is FS or HO
    FS = false;
    HO = false;
    for mii = 1:size(model_types,1)
      if strcmpi(model_types{mii,1},m)
        if strcmpi(model_types{mii,2},'FS')
          FS = true;
        else
          HO = true;
        end
      end
    end
    if ~(FS || HO)
      for mii = 1:size(model_types,1)
        if strcmpi(model_types{mii,1}(1:3),m(1:3))
          if strcmpi(model_types{mii,2},'FS')
            FS = true;
          else
            HO = true;
          end
        end
      end
    end
    if ~(FS || HO)
      % Unknown model?
      continue
    end

    % Add to data ranges for HO/FS models
    up = interp1(x,u,x_FS);
    if FS
      u_FS(:,end+1) = up;
    else
      u_HO(:,end+1) = up;
    end
    
  end
  
  % Missing data points
  m = true(size(x_FS));
  for i = 1:length(x_FS)
    if sum(isnan(u_FS(i,:)))+sum(isnan(u_HO(i,:)))>0
      m(i) = false;
    end
  end
  u_FS = u_FS(m,:);
  u_HO = u_HO(m,:);
  x_FS = x_FS(m);
  
  % ISMIP-HOM ensemble data
  uav_FS = mean(u_FS,2);
  uav_HO = mean(u_HO,2);
  sigma_FS = zeros(size(x_FS));
  sigma_HO = zeros(size(x_FS));
  for i = 1:size(u_FS,1)
    sigma_FS(i) = std(u_FS(i,:));
    sigma_HO(i) = std(u_HO(i,:));
    if isnan(sigma_FS(i)); sigma_FS(i) = 0; end
    if isnan(sigma_HO(i)); sigma_HO(i) = 0; end
  end
  umin_FS = uav_FS - sigma_FS;
  umax_FS = uav_FS + sigma_FS;
  umin_HO = uav_HO - sigma_HO;
  umax_HO = uav_HO + sigma_HO;
  
  xdata = [x_FS;flipud(x_FS)];
  ydata = [umin_FS;flipud(umax_FS)];
  set(patch_FS,'xdata',xdata,'ydata',ydata)
  
  xdata = [x_FS;flipud(x_FS)];
  ydata = [umin_HO;flipud(umax_HO)];
  set(patch_HO,'xdata',xdata,'ydata',ydata)
  
  set(line_HO,'xdata',x_FS,'ydata',uav_HO)
  set(line_FS,'xdata',x_FS,'ydata',uav_FS)
  
end

%% Calculate transect matrix for UFEMISM mesh

filename = [foldername '/test_ISMIP_HOM_A_160km_output.nc'];
mesh = read_mesh_from_file( filename);

p = [mesh.xmin,-mesh.ymax/4];
q = [mesh.xmax,-mesh.ymax/4];
n = 100;
pp = [linspace( p(1),q(1),n)', linspace( p(2),q(2),n)'];

M_trans = sparse( n,mesh.nTri);

vi = 1;
for i = 1: n
  
  vi = find_containing_vertex( mesh, pp( i,:), vi);
  
  x = pp( i,1);
  y = pp( i,2);
  
  n_c = mesh.niTri( vi);
  x_c = zeros( n_c,1);
  y_c = zeros( n_c,1);
  
  for iti = 1: mesh.niTri( vi)
    ti = mesh.iTri( vi,iti);
    x_c( iti) = mesh.TriGC( ti,1);
    y_c( iti) = mesh.TriGC( ti,2);
  end
  
  [Nf_c, Nfx_c, Nfy_c] = calc_shape_functions_2D_stag_1st_order( x, y, n_c, x_c, y_c);
  
  for iti = 1: mesh.niTri( vi)
    ti = mesh.iTri( vi,iti);
    M_trans( i,ti) = Nf_c( iti);
  end
  
end

%% Plot UFEMISM results

for xi = 1: 6
  
  L = Ls( xi);

  filename = [foldername '/test_ISMIP_HOM_A_' num2str( L) 'km_output.nc'];

  % Read data
  u_3D_b_SIASSA = ncread( filename,'u_3D_b_SIASSA');
  u_3D_b_DIVA   = ncread( filename,'u_3D_b_DIVA');
  u_3D_b_BPA    = ncread( filename,'u_3D_b_BPA');
  
  u_surf_SIASSA = u_3D_b_SIASSA( :,1);
  u_surf_DIVA   = u_3D_b_DIVA(   :,1);
  u_surf_BPA    = u_3D_b_BPA(    :,1);
  
  % Calculate transects
  u_trans_SIASSA = M_trans * u_surf_SIASSA;
  u_trans_DIVA   = M_trans * u_surf_DIVA;
  u_trans_BPA    = M_trans * u_surf_BPA;
  
  % Plot
  if (xi == 1)
    i = 1;
    j = 1;
  elseif (xi == 2)
    i = 1;
    j = 2;
  elseif (xi == 3)
    i = 1;
    j = 3;
  elseif (xi == 4)
    i = 2;
    j = 1;
  elseif (xi == 5)
    i = 2;
    j = 2;
  elseif (xi == 6)
    i = 2;
    j = 3;
  end
  
  line('parent',H2.Ax{ i,j},'xdata',linspace(-0.5,1.5,n),'ydata',u_trans_SIASSA,'color',c_SIASSA,'linewidth',3)
  line('parent',H2.Ax{ i,j},'xdata',linspace(-0.5,1.5,n),'ydata',u_trans_DIVA  ,'color',c_DIVA  ,'linewidth',3)
  line('parent',H2.Ax{ i,j},'xdata',linspace(-0.5,1.5,n),'ydata',u_trans_BPA   ,'color',c_BPA   ,'linewidth',3)

end

legend( H2.Ax{ 1,1},'Higher-Order','HO mean','Full-Stokes','FS mean','SIA/SSA','DIVA','BPA','location','northwest')

function plot_mesh_data_b( ax, mesh, d)
  patch('parent',ax,'vertices',mesh.V(1:mesh.nV,:),'faces',mesh.Tri(1:mesh.nTri,:),...
    'facecolor','flat','facevertexcdata',d,'edgecolor','none');
  set(ax,'xlim',[mesh.xmin,mesh.xmax],'ylim',[mesh.ymin,mesh.ymax]);
end

function [R,model_types] = read_Pattyn2008_ensemble( foldername)

model_types = {...
  'aas1','FS';...
  'aas2','FS';...
  'ahu1','HO';...
  'ahu2','HO';...
  'bds1','HO';...
  'cma1','FS';...
  'cma2','HO';...
  'dpo1','HO';...
  'fpa1','HO';...
  'fpa2','FS';...
  'fsa1','HO';...
  'ghg1','FS';...
  'jvj1','FS';...
  'lpe1','HO';...
  'mbr1','HO';...
  'mmr1','FS';...
  'mtk1','HO';...
  'oga1','FS';...
  'oso1','HO';...
  'rhi1','FS';...
  'rhi2','HO';...
  'rhi3','FS';...
  'rhi4','HO';...
  'rhi5','HO';...
  'spr1','FS';...
  'ssu1','FS';...
  'tpa1','HO';...
  'yko1','FS'};

models = dir(foldername);
models = models(3:end);

R.a160 = [];
R.a080 = [];
R.a040 = [];
R.a020 = [];
R.a010 = [];
R.a005 = [];

for mi = 1: length(models)
  
  modeldata = dir([foldername models(mi).name]);
  modeldata = modeldata(3:end);
  
  % Go over all experiments, check if this model has them.
  flds = fields(R);
  for xi = 1:length(flds)
    ex = flds{xi};
    
    for di = 1:length(modeldata)
      mdname = modeldata(di).name;
      str = [ex '.txt'];
      if length(mdname) >= length(str)
        if strcmpi(mdname(end-length(str)+1:end),str)
          % This is the experiment from this model
          
          disp(['Reading data from model ' models(mi).name ', experiment ' ex])
          
          fid = fopen([foldername models(mi).name '/' mdname]);
          temp = textscan(fid,'%s','delimiter','\n','MultipleDelimsAsOne',1); temp = temp{1};
          fclose(fid);
          
          n = length(temp);
          nx = sqrt(n);
          if nx-floor(nx)>0
            error('whaa!')
          end
          x_vec = zeros(n,1);
          y_vec = zeros(n,1);
          u_vec = zeros(n,1);
          v_vec = zeros(n,1);
          w_vec = zeros(n,1);
          for i = 1:n
            temp2 = textscan(temp{i},'%f %f %f %f %f %f %f %f');
            x_vec(i) = temp2{1};
            y_vec(i) = temp2{2};
            u_vec(i) = temp2{3};
            v_vec(i) = temp2{4};
            w_vec(i) = temp2{5};
          end
          
          R.(ex).(models(mi).name).x = reshape(x_vec,[nx,nx]);
          R.(ex).(models(mi).name).y = reshape(y_vec,[nx,nx]);
          R.(ex).(models(mi).name).u = reshape(u_vec,[nx,nx]);
          R.(ex).(models(mi).name).v = reshape(v_vec,[nx,nx]);
          R.(ex).(models(mi).name).w = reshape(w_vec,[nx,nx]);
          
          if (R.(ex).(models(mi).name).x(1,1) == R.(ex).(models(mi).name).x(end,1))
            R.(ex).(models(mi).name).x = R.(ex).(models(mi).name).x';
            R.(ex).(models(mi).name).y = R.(ex).(models(mi).name).y';
            R.(ex).(models(mi).name).u = R.(ex).(models(mi).name).u';
            R.(ex).(models(mi).name).v = R.(ex).(models(mi).name).v';
            R.(ex).(models(mi).name).w = R.(ex).(models(mi).name).w';
          end
          
        end
      end
    end
  end
  
end

end

function vi = find_containing_vertex( mesh, p, vi)
% Find the vertex whose Voronoi cell contains the point p, using a "linear search"
% Start at initial guess vi. Check all neighbours of vi, find the one
% closest to p, select that one as the new vi. Repeat until all neighbours
% of vi are further away from p than vi itself.

ncycle  = 0;
vi_prev = vi;
while (ncycle < mesh.nV)

  d = norm( mesh.V( vi,:) - p);

  dcmin = d + 10.0;
  vcmin = 0;
  for ci = 1: mesh.nC( vi)
    vc = mesh.C( vi,ci);
    if (vc == vi_prev); continue; end % This is the neighbour we just came from
    dc = norm( mesh.V( vc,:) - p);
    if (dc < dcmin)
      dcmin = dc;
      vcmin = vc;
    end
  end

  if (dcmin < d)
    vi_prev = vi;
    vi = vcmin;
  else
    return
  end

end % DO WHILE (ncycle < mesh.nV)

% If we reach this point, we didnt find the containing vertex - should not be possible, so throw an error
error('couldnt find closest vertex!')

end

function [Nf_c, Nfx_c, Nfy_c] = calc_shape_functions_2D_stag_1st_order( x, y, n_c, x_c, y_c)
% Calculate shape functions...
% ...in two dimensions...
% ...on the staggered grid (i.e. f is not known)...
% ...to 1st-order accuracy.
%
% Based on the least-squares approach from Syrakos et al. (2017).

% IMPLICIT NONE
% 
% % In/output variables:
% REAL(dp),                            INTENT(IN)    :: x, y       % The location where we want to know the gradients
% INTEGER,                             INTENT(IN)    :: n_max      % The maximum number of surrounding points
% INTEGER,                             INTENT(IN)    :: n_c        % The number  of     surrounding points where we know f
% REAL(dp), DIMENSION(n_max),          INTENT(IN)    :: x_c, y_c   % Coordinates of the surrounding points where we know f
% REAL(dp), DIMENSION(n_max),          INTENT(OUT)   :: Nf_c       % map    shape functions for the surrounding points
% REAL(dp), DIMENSION(n_max),          INTENT(OUT)   :: Nfx_c      % d/dx   shape functions for the surrounding points
% REAL(dp), DIMENSION(n_max),          INTENT(OUT)   :: Nfy_c      % d/dy   shape functions for the surrounding points
% LOGICAL,                             INTENT(OUT)   :: succeeded  % Whether or not we succeeded (if not, we need more neighbours)
% 
% % Local variables:
q = 1.5;
% REAL(dp), PARAMETER                                :: q = 1.5_dp
% INTEGER                                            :: ci
% REAL(dp), DIMENSION(n_c)                           :: dx, dy, w
% REAL(dp), DIMENSION(3,3)                           :: ATWTWA
% REAL(dp)                                           :: detATWTWA
% REAL(dp), DIMENSION(3,3)                           :: M

% Safety
if (n_c < 3); error('calc_shape_functions_2D_stag_1st_order needs at least 3 neighbours!'); end

% Calculate distances relative to [x,y]
dx = zeros( size( x_c));
dy = zeros( size( y_c));
for ci = 1: n_c
  dx( ci) = x_c( ci) - x;
  dy( ci) = y_c( ci) - y;
end

% Calculate the weights w
w = zeros( size( x_c));
for ci = 1: n_c
  w( ci) = 1.0 / (norm( [dx( ci), dy( ci)])^q);
end

% The matrix ATWTWA that needs to be inverted
ATWTWA = zeros( 3,3);
for ci = 1: n_c
  ATWTWA( 1,1) = ATWTWA( 1,1) + (w( ci)^2 * 1.0     * 1.0    );
  ATWTWA( 1,2) = ATWTWA( 1,2) + (w( ci)^2 * 1.0     * dx( ci));
  ATWTWA( 1,3) = ATWTWA( 1,3) + (w( ci)^2 * 1.0     * dy( ci));

  ATWTWA( 2,1) = ATWTWA( 2,1) + (w( ci)^2 * dx( ci) * 1.0    );
  ATWTWA( 2,2) = ATWTWA( 2,2) + (w( ci)^2 * dx( ci) * dx( ci));
  ATWTWA( 2,3) = ATWTWA( 2,3) + (w( ci)^2 * dx( ci) * dy( ci));

  ATWTWA( 3,1) = ATWTWA( 3,1) + (w( ci)^2 * dy( ci) * 1.0    );
  ATWTWA( 3,2) = ATWTWA( 3,2) + (w( ci)^2 * dy( ci) * dx( ci));
  ATWTWA( 3,3) = ATWTWA( 3,3) + (w( ci)^2 * dy( ci) * dy( ci));
end

% % Check if this matrix is singular
% detATWTWA = det( ATWTWA);
% if (abs( detATWTWA) < 1e-12)
%   % ATWTWA is singular
%   error('whaa!')
% end

% Invert ATWTWA to find M
M = inv( ATWTWA);

% Calculate shape functions
Nf_c    = zeros( size( x_c));
Nfx_c   = zeros( size( x_c));
Nfy_c   = zeros( size( x_c));
for ci = 1: n_c
  Nf_c(  ci) = w( ci)^2 * ( ...
    (M( 1,1) * 1.0    ) + ...
    (M( 1,2) * dx( ci)) + ...
    (M( 1,3) * dy( ci)));
  Nfx_c(  ci) = w( ci)^2 * ( ...
    (M( 2,1) * 1.0    ) + ...
    (M( 2,2) * dx( ci)) + ...
    (M( 2,3) * dy( ci)));
  Nfy_c(  ci) = w( ci)^2 * ( ...
    (M( 3,1) * 1.0    ) + ...
    (M( 3,2) * dx( ci)) + ...
    (M( 3,3) * dy( ci)));
end

end