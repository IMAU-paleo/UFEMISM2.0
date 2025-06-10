clc
clear all
close all

% Create NetCDF files containing the bed roughness and bed topography
% (including perturbed versions) for experiment II (based on the MISMIP+
% geometry, but with a spatially variable bed roughness).

%% Define parameters

% The resolutions at which we want to have the data
resolutions = [5,2.5] * 1e3;

% Domain size
xmin        = -400*1e3;     % x-coordinate of western  domain border                    [m]
xmax        =  400*1e3;     % x-coordinate of eastern  domain border                    [m]
ymin        =  -40*1e3;     % y-coordinate of southern domain border                    [m]
ymax        =   40*1e3;     % y-coordinate of northern domain border                    [m]

% Bed roughness parameters
phi_min     = 1.0;          % Till friction angle in the centre of the ice stream       [degrees]
phi_max     = 5.0;          % Till friction angle outside of the ice stream             [degrees]
x_c         = -50;          % x-coordinate of ice-stream centre                         [m]
y_c         = 0;            % y-coordinate of ice-stream centre                         [m]
sigma_x     = 150*1e3;      % x-direction ice-stream half-width                         [m]
sigma_y     =  15*1e3;      % y-direction ice-stream half-width                         [m]

% File names
filename_bed_roughness  = 'exp_II_bed_roughness';
filename_bed_topography = 'exp_II_topography';

for ri = 1:length( resolutions)
  
  %% Define the grid
  
  % The grid resolution
  resolution = resolutions( ri);
  
  % The resolution filename extension
  if resolution == 5000
    str_res = '_5km';
  elseif resolution == 2500
    str_res = '_2p5km';
  else
    error('undefined resolution!')
  end
  
  % Create the grid
  grid = create_grid( xmin, xmax, ymin, ymax, resolution);
  
  %% Calculate data
  
  % Bed roughness
  phi = zeros( grid.nx, grid.ny);
  
  for i = 1: grid.nx
    for j = 1: grid.ny
      phi( i,j) = calc_bed_roughness( phi_min, phi_max, x_c, y_c, sigma_x, sigma_y, grid.x( i), grid.y( j));
    end
  end
  
  % Start with uniform 100 m ice thickness (this really speeds things up,
  % otherwise the first 300 years will take a long time since the velocity
  % solver will be very slow)
  H = zeros( grid.nx, grid.ny) + 100;
  
  % Bed topography from MISMIP+
  b = calc_bed_topography_MISMIPplus( grid);
  
  % Surface elevation
  ice_density      =  910.0;
  seawater_density = 1028.0;
  s = H + max( -ice_density / seawater_density * H, b);
  
  %% Create and write to NetCDF files
  
  % =========================
  % ===== Bed roughness =====
  % =========================
  
  filename = [filename_bed_roughness str_res '.nc'];
  
  % Delete existing file
  if exist( filename,'file')
    delete( filename)
  end
  
  % NetCDF template
  f = create_NetCDF_template_bed_roughness( grid, filename);
  
  % Create file
  ncwriteschema( filename, f);
  
  % Write data
  ncwrite( filename,'x'       ,grid.x);
  ncwrite( filename,'y'       ,grid.y);
  ncwrite( filename,'phi_fric',phi   );
  
  % ==========================
  % ===== Bed topography =====
  % ==========================
  
  filename = [filename_bed_topography str_res '.nc'];
  
  % Delete existing file
  if exist( filename,'file')
    delete( filename)
  end
  
  % NetCDF template
  f = create_NetCDF_template_bed_topography( grid, filename);
  
  % Create file
  ncwriteschema( filename, f);
  
  % Write data
  ncwrite( filename,'x'       ,grid.x);
  ncwrite( filename,'y'       ,grid.y);
  ncwrite( filename,'Hi'      ,H);
  ncwrite( filename,'Hb'      ,b);
  ncwrite( filename,'Hs'      ,s);
  
end

function grid = create_grid( xmin, xmax, ymin, ymax, dx)
  % Create a square grid
  %
  % Code copied from IMAU-ICE
  
  grid.dx = dx;
      
  % Determine the number of grid cells that can fit in this domain
  xmid = (xmax + xmin) / 2;
  ymid = (ymax + ymin) / 2;
  nsx  = floor( (xmax - xmid) / grid.dx);
  nsy  = floor( (ymax - ymid) / grid.dx);

  grid.nx = 1 + 2*nsx;
  grid.ny = 1 + 2*nsy;
  
  % Fill in x and y
  grid.x = zeros( grid.nx,1);
  grid.y = zeros( grid.ny,1);
  
  for i = 1: grid.nx
    grid.x( i) = -nsx*grid.dx + (i-1)*grid.dx;
  end
  for j = 1: grid.ny
    grid.y( j) = -nsy*grid.dx + (j-1)*grid.dx;
  end
end
function phi  = calc_bed_roughness( phi_min, phi_max, x_c, y_c, sigma_x, sigma_y, x, y)
  % Calculate the spatially variable bed roughness field
  %
  % Basically a uniform value of phi_max, with a single ice stream at [x_c,y_c]
  % where it decreases to phi_min at the centre, with a half-width of
  % sigma_x, sigma_y in the respective x- and y-directions.
  
  phi = phi_max - (phi_max - phi_min) * exp( -0.5 * (((x - x_c) / sigma_x)^2 + ((y - y_c) / sigma_y)^2));
  
end
function b    = calc_bed_topography_MISMIPplus( grid)
  % Calculate the bed topography of the MISMIP+ experiment (Asay-Davis et al., 2016)
  
  b = zeros( grid.nx, grid.ny);
  
  B0     = -150;
  B2     = -728.8;
  B4     = 343.91;
  B6     = -50.57;
  xbar   = 300000;
  fc     = 4000;
  dc     = 500;
  wc     = 24000;
  zbdeep = -720;
  
  for i = 1: grid.nx
    for j = 1: grid.ny
      
      % Make sure everything is properly centred
      xp = grid.x( i) + 400000;
      yp = -40000 +  80000 * (j-1) / (grid.ny-1);
      
      % Asay-Davis et al. (2016), Eqs. 1-4
      xtilde = xp / xbar;
      Bx = B0 + (B2 * xtilde^2) + (B4 * xtilde^4) + (B6 * xtilde^6);
      By = (dc / (1 + exp(-2*(yp - wc)/fc))) + (dc / (1 + exp( 2*(yp + wc)/fc)));
      b( i,j) = max( Bx + By, zbdeep);
      
    end
  end
  
end
function f    = create_NetCDF_template_bed_roughness(  grid, filename)
  % Create a template for the bed roughness NetCDF file
  
  % Metadata
  f.Filename   = filename;
  f.Name       = '/';
  
  % Dimensions
  f.Dimensions(1).Name      = 'x';
  f.Dimensions(1).Length    = grid.nx;
  f.Dimensions(1).Unlimited = false;
  
  f.Dimensions(2).Name      = 'y';
  f.Dimensions(2).Length    = grid.ny;
  f.Dimensions(2).Unlimited = false;
  
  % Dimension variables
  
  % x
  f.Variables(1).Name         = 'x';
  f.Variables(1).Dimensions   = f.Dimensions(1);
  f.Variables(1).Size         = grid.nx;
  f.Variables(1).Datatype     = 'double';
  f.Variables(1).Attributes(1).Name  = 'long_name';
  f.Variables(1).Attributes(1).Value = 'X-coordinate';
  f.Variables(1).Attributes(2).Name  = 'units';
  f.Variables(1).Attributes(2).Value = 'm';
  f.Variables(1).ChunkSize    = [];
  f.Variables(1).FillValue    = [];
  f.Variables(1).DeflateLevel = [];
  f.Variables(1).Shuffle      = false;
  
  % y
  f.Variables(2).Name         = 'y';
  f.Variables(2).Dimensions   = f.Dimensions(2);
  f.Variables(2).Size         = grid.ny;
  f.Variables(2).Datatype     = 'double';
  f.Variables(2).Attributes(1).Name  = 'long_name';
  f.Variables(2).Attributes(1).Value = 'Y-coordinate';
  f.Variables(2).Attributes(2).Name  = 'units';
  f.Variables(2).Attributes(2).Value = 'm';
  f.Variables(2).ChunkSize    = [];
  f.Variables(2).FillValue    = [];
  f.Variables(2).DeflateLevel = [];
  f.Variables(2).Shuffle      = false;
  
  % Bed roughness variable
  
  f.Variables(3).Name         = 'phi_fric';
  f.Variables(3).Dimensions   = f.Dimensions;
  f.Variables(3).Size         = [grid.nx, grid.ny];
  f.Variables(3).Datatype     = 'double';
  f.Variables(3).Attributes(1).Name  = 'long_name';
  f.Variables(3).Attributes(1).Value = 'Till friction angle';
  f.Variables(3).Attributes(2).Name  = 'units';
  f.Variables(3).Attributes(2).Value = 'degrees';
  f.Variables(3).ChunkSize    = [];
  f.Variables(3).FillValue    = [];
  f.Variables(3).DeflateLevel = [];
  f.Variables(3).Shuffle      = false;
  
  % Final metadata
  f.Attributes = [];
  f.Groups     = [];
  f.Format     = 'classic';
  
end
function f    = create_NetCDF_template_bed_topography( grid, filename)
  % Create a template for the bed topography NetCDF file
  
  % Metadata
  f.Filename   = filename;
  f.Name       = '/';
  
  % Dimensions
  f.Dimensions(1).Name      = 'x';
  f.Dimensions(1).Length    = grid.nx;
  f.Dimensions(1).Unlimited = false;
  
  f.Dimensions(2).Name      = 'y';
  f.Dimensions(2).Length    = grid.ny;
  f.Dimensions(2).Unlimited = false;
  
  % Dimension variables
  
  % x
  f.Variables(1).Name         = 'x';
  f.Variables(1).Dimensions   = f.Dimensions(1);
  f.Variables(1).Size         = grid.nx;
  f.Variables(1).Datatype     = 'double';
  f.Variables(1).Attributes(1).Name  = 'long_name';
  f.Variables(1).Attributes(1).Value = 'X-coordinate';
  f.Variables(1).Attributes(2).Name  = 'units';
  f.Variables(1).Attributes(2).Value = 'm';
  f.Variables(1).ChunkSize    = [];
  f.Variables(1).FillValue    = [];
  f.Variables(1).DeflateLevel = [];
  f.Variables(1).Shuffle      = false;
  
  % y
  f.Variables(2).Name         = 'y';
  f.Variables(2).Dimensions   = f.Dimensions(2);
  f.Variables(2).Size         = grid.ny;
  f.Variables(2).Datatype     = 'double';
  f.Variables(2).Attributes(1).Name  = 'long_name';
  f.Variables(2).Attributes(1).Value = 'Y-coordinate';
  f.Variables(2).Attributes(2).Name  = 'units';
  f.Variables(2).Attributes(2).Value = 'm';
  f.Variables(2).ChunkSize    = [];
  f.Variables(2).FillValue    = [];
  f.Variables(2).DeflateLevel = [];
  f.Variables(2).Shuffle      = false;
  
  % Ice thickness variable
  
  f.Variables(3).Name         = 'Hi';
  f.Variables(3).Dimensions   = f.Dimensions;
  f.Variables(3).Size         = [grid.nx, grid.ny];
  f.Variables(3).Datatype     = 'double';
  f.Variables(3).Attributes(1).Name  = 'long_name';
  f.Variables(3).Attributes(1).Value = 'Ice thickness';
  f.Variables(3).Attributes(2).Name  = 'units';
  f.Variables(3).Attributes(2).Value = 'm';
  f.Variables(3).ChunkSize    = [];
  f.Variables(3).FillValue    = [];
  f.Variables(3).DeflateLevel = [];
  f.Variables(3).Shuffle      = false;
  
  % Bed topography variable
  
  f.Variables(4).Name         = 'Hb';
  f.Variables(4).Dimensions   = f.Dimensions;
  f.Variables(4).Size         = [grid.nx, grid.ny];
  f.Variables(4).Datatype     = 'double';
  f.Variables(4).Attributes(1).Name  = 'long_name';
  f.Variables(4).Attributes(1).Value = 'Bedrock elevation';
  f.Variables(4).Attributes(2).Name  = 'units';
  f.Variables(4).Attributes(2).Value = 'm';
  f.Variables(4).ChunkSize    = [];
  f.Variables(4).FillValue    = [];
  f.Variables(4).DeflateLevel = [];
  f.Variables(4).Shuffle      = false;
  
  % Surface elevation variable
  
  f.Variables(5).Name         = 'Hs';
  f.Variables(5).Dimensions   = f.Dimensions;
  f.Variables(5).Size         = [grid.nx, grid.ny];
  f.Variables(5).Datatype     = 'double';
  f.Variables(5).Attributes(1).Name  = 'long_name';
  f.Variables(5).Attributes(1).Value = 'Surface elevation';
  f.Variables(5).Attributes(2).Name  = 'units';
  f.Variables(5).Attributes(2).Value = 'm';
  f.Variables(5).ChunkSize    = [];
  f.Variables(5).FillValue    = [];
  f.Variables(5).DeflateLevel = [];
  f.Variables(5).Shuffle      = false;
  
  % Final metadata
  f.Attributes = [];
  f.Groups     = [];
  f.Format     = 'classic';
  
end