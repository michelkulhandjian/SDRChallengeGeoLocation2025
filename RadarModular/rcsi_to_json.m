function rcsi_to_json(csi_mat, Nfft,f_name, bin_in, ants_used)
%%%
% Gives out a json file of RCSIs for all subcarriers and all antennas.
% Each blocks holds as many entries as the antennas and there are Nfft 
% such blockd. In the end there is a block based on the metadata showing
% the antennas used for the Radar waeform that was transmitted.
%%%

n_ant = size(csi_mat,1);
f_id = fopen(f_name,'w');
init_chars ='[{ "rCSI": ';
close_crl_brkt ='},';
fprintf(f_id, '%s',init_chars );
fprintf(f_id, '%s\t\n', "[");
Nsc = Nfft;
% Loop over subcarriers
for isc=1:Nsc
    fprintf(f_id, '\t\t\t%s\n', "[");
    % Loop over antennas
    for iant = 1:n_ant
        re_im = [real(csi_mat(iant, isc)) ...
            imag(csi_mat(iant,isc))];
        jenced = jsonencode(re_im);
        if iant < n_ant
            fprintf(f_id, '\t\t\t\t%s,\n', jenced);
        else
            fprintf(f_id, '\t\t\t\t%s\n', jenced);
        end
    end
        if isc < Nsc
            fprintf(f_id, '\t\t\t%s\n', "],");
        else
            fprintf(f_id, '\t\t\t%s\n', "]");
        end
end

fprintf(f_id, '\t%s\n', "]");
fprintf(f_id, '%s\n',close_crl_brkt );
% Add used Radar waveform and antennas used
fprintf(f_id, '\n{\"%s\": %s }', bin_in, jsonencode(ants_used));
fprintf(f_id, '\n%s\n', "]");
fclose(f_id);
end