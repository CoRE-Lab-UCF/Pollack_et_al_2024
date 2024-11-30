close all
clear all
clc

%% SCRIPT TO CREATE THE BOUNDARY CONDITIONS FILES FOR SFINCS FROM THE EXTREME VALUE ANALYSIS

% To be changed to the directory where the outputs from the Extreme Value
% Analysis are saved
cd('**PATH DATA**')

folder_out_rain='**PATH WHERE THE RF SFINCS FILES WILL BE SAVED**';
folder_out_wl='**PATH WHERE THE WL SFINCS FILES WILL BE SAVED**';

data=importdata('Design_events_Cropped_Irene_.mat');
% data=importdata('Design_events_Cropped_Irene.mat');
coor=importdata('Lat_lon_Time_Irene.mat');

v2struct(coor); clear coor

lat=sort(latt,"descend");

[LON,LAT]=meshgrid(lon,lat);

% Converting Lat Lon to X Y NAD83 UTM
[A,R]=readgeoraster('**PATH_Data_Flood_Modeling**\Gloucester_street_light_utm.tif');
proj = R.ProjectedCRS;
proj.GeographicCRS.Name

[X,Y] = projfwd(proj,LAT(:),LON(:)); % x & y have to be row vectors!!

X=reshape(X,[length(lat),length(lon)]);
Y=reshape(Y,[length(lat),length(lon)]);

x=linspace(X(1,1),X(1,end),length(lon)); %x=x';
y=linspace(Y(1,1),Y(end,1),length(lat)); %y=y';

RPs=[1,2,5,10,15,20,25,50,75,100,200,500];

for rp=1:length(RPs)
    
    folder_out_rain_rp=strcat(folder_out_rain,'RP_',num2str(RPs(rp),'%03.f'));
    folder_out_wl_rp=strcat(folder_out_wl,'\','RP_',num2str(RPs(rp),'%03.f'));

    name_ncFile=strcat('Rain_',num2str(RPs(rp),'%03.f'),'.nc');
    name_wl_file=strcat('WL_',num2str(RPs(rp),'%03.f'),'.bzs');
    
    data(rp).Rainfall(isnan(data(rp).Rainfall))=-999;

% - ampr input matrix dimensions assumed to be (t,y,x)

for i=1:length(y)
    for j=1:length(x)
        for k=1:length(Time)
            PR_nc(k,i,j)=data(rp).Rainfall(i,j,k);
        end
    end
end



%% Times 
tref=datetime(Time(1),'Format','yyyyMMdd HHmmSS','convertFrom','datenum')

formatOut='yyyy-mm-dd HH:MM:SS';
refdate=datestr(tref,formatOut)

% TIMES IN PRECIPITATION IN MINUTES SINCE TREF!!!! (For WLs is in Seconds!)

time=0:60:60*(length(Time)-1);

EPSGcode=26918; 
UTMname = 'UTM18N';


sfincs_write_netcdf_amprfile(strcat(folder_out_rain_rp,'\', name_ncFile), x, y, EPSGcode, UTMname, refdate, time, PR_nc)


% Time
tref=datetime(Time(1),'Format','yyyyMMdd HHmmSS','convertFrom','datenum')

formatOut='yyyy-mm-dd HH:MM:SS';
refdate=datestr(tref,formatOut)


time_wl=0:3600:3600*(length(Time)-1);

wl_bc=[time_wl', data(rp).WL];

ascbzsfile=inputdlg({'Name for Boundary Condition WL file'},'Name BZS file',[1 50],{'sfincs'});

dlmwrite(strcat(folder_out_wl_rp,'\',name_wl_file),wl_bc,'delimiter',' ','precision','%8.2f');

end
