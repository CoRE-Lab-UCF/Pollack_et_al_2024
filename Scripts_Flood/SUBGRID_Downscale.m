clear variables;close all;

%% CODE MAARTEN FOR EXTRACTING WATER DEPTH SUBGRID RUNS

% Directories where data is saved -> You'll need to change this
output_folder = '**PATH WHERE SFINCS OUTPUT IS SAVED';
output_name=strcat('Max_WD_.mat')

cd(output_folder)

%% To be done if DEM is not in .mat format
dem_file_mat='dem_subgrid_1m_nbd.mat';

% Read input file and make grid
disp('Reading sfincs.inp ...');
inp=sfincs_read_input([output_folder 'sfincs.inp'],[]);
[xg,yg,xz,yz]=sfincs_make_grid(inp.x0,inp.y0,inp.dx,inp.dy,inp.mmax,inp.nmax,inp.rotation);

% Read output zsmax
disp('Reading output file ...');
zsmax=nc_varget([output_folder 'sfincs_map.nc'], 'zsmax');
zsmax=squeeze(zsmax);
% Read in minimum subgrid bed levels
zbmin=nc_varget([output_folder 'sfincs_map.nc'], 'zb');

% New code from Maarten not to mess up coordinates DEM
if inp.nmax==size(zsmax,2) && inp.mmax==size(zsmax,1)
    % Probably a slightly older SFINCS version (<2023) -> need to transpose
    zsmax=zsmax';
    zbmin=zbmin';
end

% Read DEM
disp('Reading DEM ...');
dem=load([output_folder dem_file_mat]);

% Flood depth below hmin will be set to NaN
hmin=0.1;
% Flood depth at bed level lower than zmin will be set to NaN
zmin=1.0;

% Get water levels on DEM grid
disp('Downscaling ...');
index_file='indices_nodata.dat'; % File name for index file (this will be made once in the current working directory)
high_zs = downscale_zs(xg, yg, zsmax, zbmin, dem.x, dem.y, hmin, index_file);

% Flood depth (subtract DEM from high-res water level)
h=high_zs-dem.z;

% Set flood depth <hmin to NaN
h(h<hmin)=NaN;

% Set flood depth at bed level <zmin to NaN
h(dem.z<zmin)=NaN;

[X,Y]=meshgrid(dem.x,dem.y);

%% Saving results to .mat file (can be converted to ArcGIS ascii file with other script)

save([output_folder output_name],'X','Y','h')
