function [latLon] = enuTogeo(pos, rP)

wgs84 = wgs84Ellipsoid;
lat0 = rP.geoData(1, 1);
lon0 = rP.geoData(2, 1);

x = pos(1); y = pos(2);
if length(pos) == 3
z = pos(3);
else
    z = 0;
end
		 
[lat, lon, ~] = enu2geodetic(x, y, z, lat0, lon0, 0, wgs84);

latLon = [lat, lon]';