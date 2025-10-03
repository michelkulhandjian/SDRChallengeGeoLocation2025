function plotDetRes(res)

loglogFlag    = true;
lnAntsperGrp  = length(res(1).pd);
lWindow        = length(res(1).pd{1});
lRes           = length(res);
cmap           = jet(lRes+3);

nAntsperGrp = res(1).nAntsperGrp;
WINDOW      = res(1).WINDOW;
%data        = res(1).dataset;

label = 1:lRes;
legendCell = cell(lRes,1);

for a = 1: lnAntsperGrp
    for w = 1: lWindow
        figSig = figure;
        figROC = figure;
        indC = 1;
        for dd = 1: lRes
            sigpMSs = res(dd).sigpMSs{a}{w};
            noisepMSs = res(dd).noisepMSs{a}{w};
            noisepMSs = noisepMSs(1:length(sigpMSs));
            figure(figSig); plot(sigpMSs, 'Color', cmap(indC,:), 'LineWidth', 2, 'LineStyle', '-');
            hold on;
            plot(noisepMSs, 'Color', cmap(indC,:), 'LineWidth', 2, 'LineStyle', '--');

            if loglogFlag
                figure(figROC); loglog(res(dd).pfa{a}{w}, res(dd).pd{a}{w}, 'Color', cmap(indC,:), 'LineWidth', 2);
            else
                figure(figROC); semilogx(res(dd).pfa{a}{w}, res(dd).pd{a}{w}, 'Color', cmap(indC,:), 'LineWidth', 2);
            end
            hold on;
            legendCell{indC} =  int2str(label(indC));
            indC = indC+1;
        end
        figure(figSig); title(['Detection Output, nAntsperGrp = ' int2str(nAntsperGrp(a)) ', W = ' int2str(WINDOW(w))]);
        xlabel('Samples'); ylabel('Values'); grid on;
        legend(legendCell,'Location','Best');

        figure(figROC); title(['ROC, nAntsperGrp = ' int2str(nAntsperGrp(a)) ', W = ' int2str(WINDOW(w))]);
        xlabel('Pfa'); ylabel('Pd'); grid on;
        legend(legendCell,'Location','Best');

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
end