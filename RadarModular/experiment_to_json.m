function experiment_to_json( n_ant_per_grp, n_freq_grp, act_ants, f_name)

if ~exist('f_name','var')
    f_name = 'experiment.json';
end
if ~exist('n_ant_per_grp','var')
    n_ant_per_grp = 14;
end
if ~exist('n_freq_grp','var')
    n_freq_grp = 6;
end
if ~exist('act_ants','var')
    act_ants = ['[ 2, 3, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, ' ...
        '19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 44, 45, 46, 47, 48, 49, 50, 51, 52, ' ...
        '53, 54, 55, 56, 57, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71 ]'];
end

f_id = fopen(f_name,'w');
init_brkt ='{';
end_brkt ='}';
num_ant_per_grp_str ='"num_ant_per_group": ';
num_freq_grp_str ='"num_freq_groups": ';
active_ant_str  = '"active_antennas": ';
comma           = ',';
fprintf(f_id, '%s\n',init_brkt );
fprintf(f_id, '\t%s',num_ant_per_grp_str );
fprintf(f_id, '%d',n_ant_per_grp );
fprintf(f_id, '%s\n',comma );
fprintf(f_id, '\t%s',num_freq_grp_str );
fprintf(f_id, '%d',n_freq_grp );
fprintf(f_id, '%s\n',comma );
fprintf(f_id, '\t%s',active_ant_str );
fprintf(f_id, '%s\n',act_ants );
fprintf(f_id, '%s\n',end_brkt );
fclose(f_id);
end

%{
 %   "num_ant_per_group": 14,
 %   "num_freq_groups": 6,
 %   "active_antennas": [ 2,  3,  6,  7,  8,  9,  10, 11, 12, 13, 14, 15, 16, 17, 18,
 %   19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 44, 45, 46, 47, 48, 49, 50, 51, 52,
 %   53, 54, 55, 56, 57, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71 ]
%}


