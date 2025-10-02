% Write to a file
function writeToFile (x, filename, dataType, dataScaling)

% Make sure that you pass raw vector (with column vector will get you wrong
% values)
write_data = zeros(length(x) * 2, 1);
write_data(1:2:end) = real(x);
write_data(2:2:end) = imag(x);
%checkmat = int16(storedInteger(fi(write_data, 1, 16, 11)));
write_data = int16(double(write_data) ./ dataScaling);
%They are not equal (very minor, rounding differences), max difference +1
%disp("Is equal? " + isequal(checkmat, write_data))
fidw = fopen(filename,'w' );
fwrite(fidw, write_data, dataType);
fclose(fidw); clear fidw
end

