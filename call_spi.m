clear
close all


%% ************************ Modify this region: ************************
% file
infile='/scratch/ault/b40.lm1850-2005.1deg.001.cam2.h0.PRECT.185001-200512.land_only.nc';
odir = '/scratch/ault/';
dashq=find(infile=='/');
if ~isempty(dashq);
    fprfx=infile(dashq(end)+1:end-3);
else
    prfx=infile(1:end-3);
end

ofilename=[odir fprfx '.spi.nc'];

%lat, lon, time (might be different for different nc files)
%e.g.:
% pvarmame='prcp';
% tvarname='T';
% lonname='X';
% latname='Y';
pvarmame='PRECT';
tvarname='time';
lonname='lon';
latname='lat';


%precipitation conversion (need mm/mo)
prscl=60*60*24*30*1000;

% Index value (NOT year) of desired start/stop
calib_start_yr=nan;
calib_end_yr=nan;

% length (in months) of desired SPIs
spilen=[1,3,6,9,12,24,48,120];

%% ************************  STOP modifying ... ***********************

% open file
ncidin=netcdf.open(infile,0);

% find correct indices for each variable
time_indx=netcdf.inqVarID(ncidin,tvarname);
lon_indx=netcdf.inqVarID(ncidin,lonname);
lat_indx=netcdf.inqVarID(ncidin,latname);
pr_indx=netcdf.inqVarID(ncidin,pvarmame);

% now import 
lat=netcdf.getVar(ncidin,lat_indx);
lon=netcdf.getVar(ncidin,lon_indx);
time=netcdf.getVar(ncidin,time_indx);
time_units=netcdf.getAtt(ncidin,time_indx,'units');
pr=netcdf.getVar(ncidin,pr_indx);


% convert units:
pr = pr*prscl;


% get rid of garbage, if present (might be different for different files)
pr(pr<0)=nan;
pr(pr>10000)=nan;

% Index years (used to make sure all year values are standardized, regardless of input)
beg_indexyr=1;
end_indexyr=length(time)/12;
index_years=beg_indexyr:end_indexyr;
true_years=time(1:12:end);

% Defaults:
if (isnan(calib_start_yr))
    calib_start_yr=1;
end

if (isnan(calib_end_yr))
    calib_end_yr=length(true_years);
end

% set up yr | mo vector fo SPI inpout file
yrs=reshape(repmat(index_years,12,1),[],1); 
yrmovec=[yrs repmat((1:12)',length(yrs)/12,1)];

% set up 4D array for SPI:
dims=size(pr);
nspi=length(spilen);
spi4d=nan([dims nspi]);

netcdf.close(ncidin)

% loop through all grids
for i =1:dims(1);
    for j= 1:dims(2);
        prcp1d = squeeze(pr(i,j,:));
        if sum(~isnan(prcp1d)) == length(pr) % all values must be present
            fid=fopen('prcptmp.txt','w');
            fprintf(fid,'%3s %3s %3s',['yr ','mo ','pr ']);
            fprintf(fid,'\n');
            fprintf(fid,'%i     %i      %i\n',reshape([yrmovec round(100*prcp1d)]',[],3));
            fclose(fid);
            eval(['!./spi  ' num2str(spilen) ' -bc ' num2str(calib_start_yr) ...
                  ' -ec ' num2str(calib_end_yr) ' < prcptmp.txt >&err.log']);

            % spi writes to ofile.txt, so
            % now we import ofile.txt, convert -99.00 to NAN and store in SPI4D
            Xout=load('ofile.txt');
            Xout(Xout<=-99.0)=nan;

            if (any(Xout==9999))
                disp('Warning: missing values in SPI output...')
                Xout(Xout==9999)=nan;
            end
            

            % Figure out index (in time) of the columns in Xout
            % Usually, this will just be the month index given in Xout(1,2), but if entire
            % chunks of years are skipped (for instance with long SPI lengths: > 12) then
            % we'll need to skip years too.
            qtstart=(find(Xout(1,1)==index_years)-1)*12+Xout(1,2);

            % loop through SPIs and store in correct place
            for k=1:nspi
                spi4d(i,j,qtstart:length(time),k)=Xout(:,2+k);
            end            
        end
    end
    % show progress
    disp([num2str(round((i-1)/dims(1)*10000)/100) '%'])
end

% change NANs to -99.9
ammsng=-99.9;
spi4d(isnan(spi4d))=ammsng;

%% netcdf write part
varname='spi4d';
latname='lat';
lonname='lon';
timename='time';

eval(['!rm ' ofilename])

mode = netcdf.getConstant('CLASSIC_MODEL');
mode = bitor(mode,netcdf.getConstant('CLASSIC_MODEL'));
ncid = netcdf.create(ofilename,mode);

lat_dimID = netcdf.defDim(ncid,'lat',length(lat));
lon_dimID = netcdf.defDim(ncid,'lon',length(lon));
spi_dimID = netcdf.defDim(ncid,'spilen',length(spilen));
time_dimID = netcdf.defDim(ncid,'time',length(time));

% SPI + attributes
var_varid = netcdf.defVar(ncid,'spi4d','double',[lon_dimID,lat_dimID,time_dimID,spi_dimID]);
netcdf.putAtt(ncid,var_varid,'missing_value',ammsng)
netcdf.putAtt(ncid,var_varid,'long_name','Standardized Precipitation Index')
netcdf.putAtt(ncid,var_varid,'units','stddev')


% Latitude + attributes
lat_varid = netcdf.defVar(ncid,'lat','double',lat_dimID);
netcdf.putAtt(ncid,lat_varid,'units','degrees_north')
netcdf.putAtt(ncid,lat_varid,'long_name','lattitude')

% Longitude + attributes
lon_varid = netcdf.defVar(ncid,'lon','double',lon_dimID);
netcdf.putAtt(ncid,lon_varid,'units','degrees_east')
netcdf.putAtt(ncid,lon_varid,'long_name','longitude')

% spilen + atts:
spi_varid = netcdf.defVar(ncid,'spilen','double',spi_dimID);
netcdf.putAtt(ncid,spi_varid,'long_name','SPI lengths');
netcdf.putAtt(ncid,spi_varid,'units','months');

% time + atts
time_varid = netcdf.defVar(ncid,'time','double',time_dimID);
netcdf.putAtt(ncid,time_varid,'long_name','time');
netcdf.putAtt(ncid,time_varid,'units',time_units);

% global attributes:
title = ['SPI computed from precipitation (' pvarmame ').'];
notes = ['SPI values computed using C++ code provided by the "Green Leaf" project (http://greenleaf.unl.edu/downloads/). Original code was written in Fortran and provided online by Guttman et al., 1999 (spinew.f). Method uses the partial gamma distribution to convert monthly precipitation into a normally distributed variable.'];
history = ['Created from precipitation values in ' infile ', with spi-lengths= [' num2str(spilen) '] months'];
Reference = 'Guttman, N. B. (1999), Accepting the standardized precipitation index: A calculation algorithm, Journal of the American Water Resources Association, 35(2), 311-322.';

netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'title',title)
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'notes',notes)
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'history',history)
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'Reference',Reference)

netcdf.endDef(ncid)

% Now actually write data to file
netcdf.putVar(ncid,var_varid,spi4d);
netcdf.putVar(ncid,lat_varid,lat);
netcdf.putVar(ncid,lon_varid,lon);
netcdf.putVar(ncid,time_varid,time);
netcdf.putVar(ncid,spi_varid,spilen);

netcdf.close(ncid)
