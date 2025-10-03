% Plot Figures
function plotROC ( res, rP )

loglogFlag    = true;
SNR           = res.SNR;            % SNRs
%SNR           = [1:6];
WINDOW        = res.ROCWindow;      %[16, 30, 60];
ANTSPERGRP    = res.ROCnAntsperGrp; % Antennas per group
NSPELEM       = res.ROCNspelem;     % Sparseness
%WINDOW        = WINDOW(1);

Nsnr          = length(SNR);
Nwin          = length(WINDOW);
NantsperGrp   = length(ANTSPERGRP);
Nsp           = length(NSPELEM);


% Plot

if rP.ROCFig  == 1
    indAntInt  = 1;
    indW       = 1;
    indNsp     = 1;
    param_str = 'SNR';            % SNRs
    %waveforms = {'BIN1-A', 'BIN1-B', 'BIN2-A', 'BIN2-B', 'BIN3-A', 'BIN3-B'};
    param_fixed_str = 'SP';       % Fixed parameter is SF (SP)  % FOR FIG TITLE
    figure;
    cmap        = hsv(Nsnr); %gray(20);   % number of sweep params
    cnt = 0;

    for indSNR =1: Nsnr
        %for indW = 1:Nwin
        snr = SNR(indSNR);
        cnt = cnt + 1;
        if loglogFlag
            loglog(res(indSNR).pfa{indAntInt}{indW}, res(indSNR).pd{indAntInt}{indW}, 'Color', cmap(cnt,:), 'LineWidth', 2, 'LineStyle', '--');
        else
            semilogx(res(indSNR).pfa{indAntInt}{indW}, res(indSNR).pd{indAntInt}{indW}, 'Color', cmap(cnt,:), 'LineWidth', 2, 'LineStyle', '--');
        end
        hold on;
        legendCell{cnt} = sprintf('%s: %0.2f ', param_str, snr);
        %legendCell{cnt} = sprintf('%s',  waveforms{indSNR});
    end

    if (1)

        yl = yline(0.9,'-.','90 %', 'HandleVisibility','off');
        yl.LabelHorizontalAlignment = 'left';
        yl.LabelVerticalAlignment = 'middle';
        x1 = xline(0.05,'-.','5 %', 'HandleVisibility','off');
        x1.LabelHorizontalAlignment = 'center';
        x1.LabelVerticalAlignment = 'middle';
        x2 = xline(0.10,'-.','10 %', 'HandleVisibility','off');
        x2.LabelHorizontalAlignment = 'center';
        x2.LabelVerticalAlignment = 'middle';
        x3 = xline(0.15,'-.','15 %', 'HandleVisibility','off');
        x3.LabelHorizontalAlignment = 'center';
        x3.LabelVerticalAlignment = 'middle';
    end

    title(sprintf('%s: %d,  Ant: %d, R.Win: %d', param_fixed_str, NSPELEM(indNsp),  rP.ROCnAntsperGrp(indAntInt), rP.ROCWindow(indW)));
    xlabel('Pfa'); ylabel('Pd'); %grid on;
    axis([(10^-4) 1 0.7 1]);
    legend(legendCell,'Location','Best');
end

% 3. Sweeping SNR, W, AntPerGrp, fix sp  - 3x3 ROCs (several SNRs each plot)
if rP.ROCFig  == 3
    indNsp  = 1;
    param_str = 'SNR';      % SNRs
    param_fixed_str = 'SP';       % Fixed parameter is SF (SP)  % FOR FIG TITLE
    figure;
    cmap        = hsv(Nsnr); %gray(20);   % number of sweep params
    cnt = 0;

    % Per antenna value
    for indAntInt = 1:NantsperGrp
        ant = ANTSPERGRP(indAntInt);
        % Per WIN
        for indW = 1:Nwin
            curW = WINDOW(indW);
            cnt = cnt + 1;
            subplot(3,3,cnt);
            % Sweep
            for indSNR = 1:Nsnr
                snr = SNR(indSNR);
                if loglogFlag
                    loglog(res(indSNR).pfa{indAntInt}{indW}, res(indSNR).pd{indAntInt}{indW}, 'Color', cmap(indSNR,:), 'LineWidth', 2, 'LineStyle', '--');
                else
                    semilogx(res(indSNR).pfa{indAntInt}{indW}, res(indSNR).pd{indAntInt}{indW}, 'Color', cmap(indSNR,:), 'LineWidth', 2, 'LineStyle', '--');
                end
                hold on;
                legendCell{indSNR} = sprintf('%s: %0.2f ', param_str, snr);

            end

            if (1)
                yl = yline(0.9,'-.','90 %', 'HandleVisibility','on');
                yl.LabelHorizontalAlignment = 'left';
                yl.LabelVerticalAlignment = 'middle';
                x1 = xline(0.05,'-.','5 %', 'HandleVisibility','on');
                x1.LabelHorizontalAlignment = 'center';
                x1.LabelVerticalAlignment = 'middle';
                x2 = xline(0.10,'-.','10 %', 'HandleVisibility','on');
                x2.LabelHorizontalAlignment = 'center';
                x2.LabelVerticalAlignment = 'middle';
                x3 = xline(0.15,'-.','15 %', 'HandleVisibility','on');
                x3.LabelHorizontalAlignment = 'center';
                x3.LabelVerticalAlignment = 'middle';
            end
            %title(sprintf(num2str([WINDOW SNR],'Win: %d, SNR: %0.1f dB')));
            title(sprintf('%s: %d, WIN: %d dB, Ant: %d', param_fixed_str, NSPELEM(indNsp), curW, ant));
            xlabel('Pfa'); ylabel('Pd'); %grid on;
            axis([(10^-4) 1 0.7 1]);
            %legend(legendCell,'Location','Best');

        end
    end
    legend(legendCell);  % Ignore the empty slots due to absence of MF (dirty hack)
    %legend(legendCell{[1:3,5,6,8,9,11,12],:});
    %legend(legendCell{[1:3,5,6],:});  % Ignore the empty slots due to absence of MF (dirty hack)
    fig = gcf;
    %fig.Position(3) = fig.Position(3) + 0*300;
    % add legend
    Lgnd = legend('show');
    Lgnd.Position(1) = 0.651;
    Lgnd.Position(2) = 0.141;

end

% 5. Sweeping W, SNR, AntPerGrp, fix sp  - 3x3 ROCs (several Ws each plot) - prev
if rP.ROCFig  == 5
    indNsp  = 1;
    param_str = 'WIN';      % WINs
    param_fixed_str = 'SP';       % Fixed parameter is SF (SP)  % FOR FIG TITLE
    figure;
    cmap        = hsv(Nwin); %gray(20);   % number of sweep params
    cnt = 0;

    % Per antenna value
    for indAntInt = 1:NantsperGrp
        ant = ANTSPERGRP(indAntInt);
        % Per SNR
        for indSNR = 1:Nsnr
            snr = SNR(indSNR);
            cnt = cnt + 1;
            subplot(3,3,cnt);
            % Sweep
            for indW = 1:Nwin
                curW = WINDOW(indW);
                %loglog(res(snrIdx).pfaMSs{antIdx}, res(snrIdx).pdMSs{antIdx}, 'Color', cmap(iswp,:), 'LineWidth', 2, 'LineStyle', '--');

                if loglogFlag
                    loglog(res(indSNR).pfa{indAntInt}{indW}, res(indSNR).pd{indAntInt}{indW}, 'Color', cmap(indW,:), 'LineWidth', 2, 'LineStyle', '--');
                else
                    semilogx(res(indSNR).pfa{indAntInt}{indW}, res(indSNR).pd{indAntInt}{indW}, 'Color', cmap(indW,:), 'LineWidth', 2, 'LineStyle', '--');
                end
                hold on;
                legendCell{indW} = sprintf('%s: %0.2f ', param_str, curW);

            end

            if (1)
                yl = yline(0.9,'-.','90 %', 'HandleVisibility','on');
                yl.LabelHorizontalAlignment = 'left';
                yl.LabelVerticalAlignment = 'middle';
                x1 = xline(0.05,'-.','5 %', 'HandleVisibility','on');
                x1.LabelHorizontalAlignment = 'center';
                x1.LabelVerticalAlignment = 'middle';
                x2 = xline(0.10,'-.','10 %', 'HandleVisibility','on');
                x2.LabelHorizontalAlignment = 'center';
                x2.LabelVerticalAlignment = 'middle';
                x3 = xline(0.15,'-.','15 %', 'HandleVisibility','on');
                x3.LabelHorizontalAlignment = 'center';
                x3.LabelVerticalAlignment = 'middle';
            end
            %title(sprintf(num2str([WINDOW SNR],'Win: %d, SNR: %0.1f dB')));
            title(sprintf('%s: %d, SNR: %d dB, Ant: %d', param_fixed_str, NSPELEM(indNsp), snr, ant));
            xlabel('Pfa'); ylabel('Pd'); %grid on;
            axis([(10^-4) 1 0.7 1]);
            %legend(legendCell,'Location','Best');

        end
    end
    legend(legendCell);  % Ignore the empty slots due to absence of MF (dirty hack)
    %legend(legendCell{[1:3,5,6,8,9,11,12],:});
    %legend(legendCell{[1:3,5,6],:});  % Ignore the empty slots due to absence of MF (dirty hack)
    fig = gcf;
    %fig.Position(3) = fig.Position(3) + 0*300;
    % add legend
    Lgnd = legend('show');
    Lgnd.Position(1) = 0.651;
    Lgnd.Position(2) = 0.141;

end

% 7. Sweeping AntPerGrp, SNR, W, fix sp  - 3x3 ROCs (several AntPerGrps each plot)
if rP.ROCFig  == 7
    indNsp  = 1;
    param_str = 'Ant';      % Ants
    param_fixed_str = 'SP';       % Fixed parameter is SF (SP)  % FOR FIG TITLE
    figure;
    cmap        = hsv(NantsperGrp); %gray(20);   % number of sweep params
    cnt = 0;


    % Per window value
    for indW = 1:Nwin
        curW = WINDOW(indW);
        % Per SNR
        for indSNR = 1:Nsnr
            snr = SNR(indSNR);
            cnt = cnt + 1;
            subplot(3,3,cnt);
            % Sweep
            for indAntInt = 1:NantsperGrp
                ant = ANTSPERGRP(indAntInt);
                if loglogFlag
                    loglog(res(indSNR).pfa{indAntInt}{indW}, res(indSNR).pd{indAntInt}{indW}, 'Color', cmap(indAntInt,:), 'LineWidth', 2, 'LineStyle', '--');
                else
                    semilogx(res(indSNR).pfa{indAntInt}{indW}, res(indSNR).pd{indAntInt}{indW}, 'Color', cmap(indAntInt,:), 'LineWidth', 2, 'LineStyle', '--');
                end
                hold on;
                legendCell{indAntInt} = sprintf('%s: %0.2f ', param_str, ant);
            end

            if (1)
                yl = yline(0.9,'-.','90 %', 'HandleVisibility','on');
                yl.LabelHorizontalAlignment = 'left';
                yl.LabelVerticalAlignment = 'middle';
                x1 = xline(0.05,'-.','5 %', 'HandleVisibility','on');
                x1.LabelHorizontalAlignment = 'center';
                x1.LabelVerticalAlignment = 'middle';
                x2 = xline(0.10,'-.','10 %', 'HandleVisibility','on');
                x2.LabelHorizontalAlignment = 'center';
                x2.LabelVerticalAlignment = 'middle';
                x3 = xline(0.15,'-.','15 %', 'HandleVisibility','on');
                x3.LabelHorizontalAlignment = 'center';
                x3.LabelVerticalAlignment = 'middle';
            end
            %title(sprintf(num2str([WINDOW SNR],'Win: %d, SNR: %0.1f dB')));
            title(sprintf('%s: %d, SNR: %d dB, WIN: %d', param_fixed_str, NSPELEM(indNsp), snr, curW));
            xlabel('Pfa'); ylabel('Pd'); %grid on;
            axis([(10^-4) 1 0.7 1]);
            %legend(legendCell,'Location','Best');

        end
    end
    legend(legendCell);  % Ignore the empty slots due to absence of MF (dirty hack)
    %legend(legendCell{[1:3,5,6,8,9,11,12],:});
    %legend(legendCell{[1:3,5,6],:});  % Ignore the empty slots due to absence of MF (dirty hack)
    fig = gcf;
    %fig.Position(3) = fig.Position(3) + 0*300;
    % add legend
    Lgnd = legend('show');
    Lgnd.Position(1) = 0.651;
    Lgnd.Position(2) = 0.141;

end

% 9. Sweeping sp, SNR, AntPerGrp, fix W - 3x3 ROCs (several sps each plot) - prev
if rP.ROCFig  == 9
    indW = 1;
    param_str = 'SP';      % SPs
    param_fixed_str = 'WIN';       % Fixed parameter is SF (SP)  % FOR FIG TITLE
    figure;
    cmap        = hsv(Nsp); %gray(20);   % number of sweep params
    cnt = 0;

    % Per antenna value
    for indAntInt = 1:NantsperGrp
        ant = ANTSPERGRP(indAntInt);
        % Per SNR
        for indSNR = 1:Nsnr
            snr = SNR(indSNR);
            cnt = cnt + 1;
            subplot(3,3,cnt);
            % Sweep
            for indNsp = 1:Nsp
                curNsp = NSPELEM(indNsp);
                %loglog(res(snrIdx).pfaMSs{antIdx}, res(snrIdx).pdMSs{antIdx}, 'Color', cmap(iswp,:), 'LineWidth', 2, 'LineStyle', '--');

                if loglogFlag
                    loglog(res(indSNR).pfa{indAntInt}{indNsp}, res(indSNR).pd{indAntInt}{indNsp}, 'Color', cmap(indNsp,:), 'LineWidth', 2, 'LineStyle', '--');
                else
                    semilogx(res(indSNR).pfa{indAntInt}{indNsp}, res(indSNR).pd{indAntInt}{indNsp}, 'Color', cmap(indNsp,:), 'LineWidth', 2, 'LineStyle', '--');
                end
                hold on;
                legendCell{indNsp} = sprintf('%s: %0.2f ', param_str, curNsp);

            end

            if (1)
                yl = yline(0.9,'-.','90 %', 'HandleVisibility','on');
                yl.LabelHorizontalAlignment = 'left';
                yl.LabelVerticalAlignment = 'middle';
                x1 = xline(0.05,'-.','5 %', 'HandleVisibility','on');
                x1.LabelHorizontalAlignment = 'center';
                x1.LabelVerticalAlignment = 'middle';
                x2 = xline(0.10,'-.','10 %', 'HandleVisibility','on');
                x2.LabelHorizontalAlignment = 'center';
                x2.LabelVerticalAlignment = 'middle';
                x3 = xline(0.15,'-.','15 %', 'HandleVisibility','on');
                x3.LabelHorizontalAlignment = 'center';
                x3.LabelVerticalAlignment = 'middle';
            end
            %title(sprintf(num2str([WINDOW SNR],'Win: %d, SNR: %0.1f dB')));
            title(sprintf('%s: %d, SNR: %d dB, Ant: %d', param_fixed_str, WINDOW(indW), snr, ant));
            xlabel('Pfa'); ylabel('Pd'); %grid on;
            axis([(10^-4) 1 0.7 1]);
            %legend(legendCell,'Location','Best');

        end
    end
    legend(legendCell);  % Ignore the empty slots due to absence of MF (dirty hack)
    %legend(legendCell{[1:3,5,6,8,9,11,12],:});
    %legend(legendCell{[1:3,5,6],:});  % Ignore the empty slots due to absence of MF (dirty hack)
    fig = gcf;
    %fig.Position(3) = fig.Position(3) + 0*300;
    % add legend
    Lgnd = legend('show');
    Lgnd.Position(1) = 0.651;
    Lgnd.Position(2) = 0.141;

end

end