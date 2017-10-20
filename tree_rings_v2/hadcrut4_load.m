%load hadcrut4 dataset

%load HADCRUT4 data
%% ************************ Modify this region: ************************
% file
infile='/Users/adam/Desktop/Treerings/data/HadCRUT.4.1.1.0.median 1.nc';
%odir = '/Users/adam/Desktop/Treerings/';
dashq=find(infile=='/');
if ~isempty(dashq);
    fprfx=infile(dashq(end)+1:end-3);
else
    prfx=infile(1:end-3);
end

%ofilename=[odir fprfx '.dat.nc'];

%lat, lon, time (might be different for different nc files)
%e.g.:
% pvarmame='prcp';
% tvarname='T';
% lonname='X';
% latname='Y';
dvarmame='temperature_anomaly';
tvarname='time';
lonname='longitude';
latname='latitude';




% Index value (NOT year) of desired start/stop
%calib_start_yr=nan;
%calib_end_yr=nan;

% length (in months) of desired SPIs
%spilen=[12,24,48,120];

%% ************************  STOP modifying ... ***********************

% open file
ncidin=netcdf.open(infile,0);

% find correct indices for each variable
time_indx=netcdf.inqVarID(ncidin,tvarname);
lon_indx=netcdf.inqVarID(ncidin,lonname);
lat_indx=netcdf.inqVarID(ncidin,latname);
d_indx=netcdf.inqVarID(ncidin,dvarmame);

% now import 
H.lat=netcdf.getVar(ncidin,lat_indx);
H.lon=netcdf.getVar(ncidin,lon_indx);
H.time=netcdf.getVar(ncidin,time_indx);
H.time_units=netcdf.getAtt(ncidin,time_indx,'units');
H.d=netcdf.getVar(ncidin,d_indx);


save('hadcrut4.mat', 'H')