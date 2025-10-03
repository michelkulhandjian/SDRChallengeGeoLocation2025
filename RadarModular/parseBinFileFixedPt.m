% Parsing the SEC data (5.11 int16)
function [output_mat] = parseBinFileFixedPt (filename)
    %Attempt to reconstruct input matrix
    % https://www.mathworks.com/help/dsp/ug/write-and-read-fixed-point-data.html
    fidr = fopen(filename, 'r');
    % input_data = fread(fidr, inf, 'int16=>int16');
    % input as int16 as double so we can scale
    input_data = fread(fidr, inf, 'int16');
    fclose(fidr); clear fidr
    fraction_length = 11;
    word_length = 16;
    fixed_pt_scaled = 2^(-fraction_length)*input_data;
    fixed_pt_scaled_complex = complex(fixed_pt_scaled(1:2:end), fixed_pt_scaled(2:2:end));
    output_mat = fi(fixed_pt_scaled_complex, 1, word_length, fraction_length);
    clear fixed_pt_scaled fixed_pt_scaled_complex fraction_length word_length dataToWrite_check
end

