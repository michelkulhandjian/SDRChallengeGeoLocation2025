% Write fixed pt int16 to file
function writeToFileFixedPt (dataToWrite, filename)
    elements = length(dataToWrite);
    % restructure as real[0], imag[0], real[1], imag[1], etc
    x_real = real(dataToWrite);
    x_imag = imag(dataToWrite);
    p(1:2:2*elements, :) = x_real;
    p(2:2:2*elements, :) = x_imag;
    % https://www.mathworks.com/help/fixedpoint/ref/embedded.fi.storedinteger.html
    store_int = storedInteger(p);
    fidw = fopen(filename,'w');
    fwrite(fidw, store_int, 'int16');
    fclose(fidw); clear fidw
    if 0
        check_mat = parseBinFileFixedPt(filename);
        disp("Match? " + isequal(dataToWrite, check_mat));
    end
end