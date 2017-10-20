%messing around with netcdf commands



ncfile = netcdf.open('/Users/adam/Desktop/Treerings/data/HadCRUT.4.1.1.0.median 1.nc',0);

[ndims, nvars, natts] = netcdf.inq(ncfile)

for i = 1:nvars
    [varname, xtype, varDimIDs, varAtts] = netcdf.inqVar(ncfile,i-1)
end
    