ncid = netcdf.open('ssi_v03r00_yearly_s1610_e2023_c20240831.nc')
A_number = netcdf.getVar(ncid,2)
netcdf.close(ncid)