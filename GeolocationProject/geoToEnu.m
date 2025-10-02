function [A] = geoToEnu (Meas, Nd, R)

if ~exist('R','var')
     R = 1;
end
%Meas = Meas';

[L, Nr] = size(Meas);

if 1

    wgs84 = wgs84Ellipsoid;

    lat0 = Meas(1, 1);
    lon0 = Meas(2, 1);
    %t0 = measurements(1, 3);
    for i = 1:Nr
        lat = Meas(1, i);
        lon = Meas(2, i);
        %t = measurements(i, 3);
        [e, n, ~] = geodetic2enu(lat, lon, 0, lat0, lon0, 0, wgs84);
        A(:, i) = [e, n]';
    end

else

    %R = 1; lat = [56 34]; lon = [-140 23];
    for i = 1:Nr
        lat = Meas(1, i);
        lon = Meas(2, i);
        %t = measurements(i, 3);
        x = R * cos(lat) .* cos(lon);
        y = R * cos(lat) .* sin(lon);
        if Nd <= 2
            A(:, i) = [x, y]';
        else
            z = R * sin(lat);
            A(:, i) = [x, y, z]';

        end

    end

end

