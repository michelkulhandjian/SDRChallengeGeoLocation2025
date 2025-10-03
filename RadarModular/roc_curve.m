% ROC Analysis
function [sim_pd, sim_pfa, thresh] = roc_curve ( h1, h0)

if isempty(h1) | isempty(h0)
    sim_pd = 0;
    sim_pfa = 0;
    return;
end
select_histogram_plot = 1;
select_roc_plot = 0;

L = length(h1);
h1a = abs(h1);
h0a = abs(h0);

thresh_low = min([h1a,h0a]);
thresh_hi  = max([h1a,h0a]);

nbins = 1000;
%fprintf( 'thresh_low %d thresh_hi %d sec \n', thresh_low, thresh_hi);
thresh_steps = linspace(thresh_low,thresh_hi,nbins);
sim_pd = zeros(1,nbins);
sim_pfa = zeros(1,nbins);
for k = 1:nbins
    thresh = thresh_steps(k);
    sim_pd(k) = sum(h1a >= thresh);
    sim_pfa(k) = sum(h0a >= thresh);
end

ind_pd1 = find(sim_pd<L)-1;
ind_pd09 = find(sim_pd<L*0.9);
ind_pfa = find(sim_pfa);
ind_pfa005 = find(sim_pfa<L*0.05);

%ind_pd1 = find(fliplr(sim_pd)>=L);
%ind_pd09 = find(fliplr(sim_pd)>L*0.9);
%ind_pfa = find(fliplr(sim_pfa));

sim_pd = sim_pd/L;
sim_pfa = sim_pfa/L;


if ind_pfa(end) < ind_pd1(1)
    ind = fix((ind_pfa(end) + ind_pd1(1))/2);
elseif ind_pfa005(1) < ind_pd1(1)
    ind = ind_pd1(1);
elseif ind_pfa005(1) < ind_pd09(1)
    ind = fix( (ind_pfa005(1)+ ind_pd09(1))/2);
else
    ind = ind_pd09(1);
end
thresh  = thresh_steps (ind );

sp = thresh_steps(2) - thresh_steps(1);
indA = find(thresh_steps >=8.5-sp & thresh_steps>=8.5-sp);
indB = find(thresh_steps >=10.5-sp & thresh_steps>=10.5-sp);
[thresh_steps(indA(1)), sim_pd(indA(1)), sim_pfa(indA(1))*100,thresh_steps(indB(1)), sim_pd(indB(1)), sim_pfa(indB(1))*100]

pfa_diff = diff(sim_pfa);
idx = (pfa_diff == 0);
sim_pfa(idx) = [];
sim_pd(idx) = [];

minpfa = 1e-6;
N = sum(sim_pfa >= minpfa);
sim_pfa = fliplr(sim_pfa(1:N)).';
sim_pd = fliplr(sim_pd(1:N)).';

nbins = 100;
binedges = linspace(thresh_low,thresh_hi,nbins);
if select_histogram_plot == 1
    figure
    histogram(h0a,binedges)
    hold on
    histogram(h1a,binedges)
    xline(thresh)
    hold off
    title('Target-Absent Vs Target-Present Histograms')
    legend('Target Absent','Target Present')
end

if select_roc_plot  == 1
    %semilogx(sim_pfa,sim_pd,'r.')
    figure
    semilogx(sim_pfa,sim_pd,'r')
    hold on
    title(' ROC Curves')
    xlabel('Pfa')
    ylabel('Pd')
    grid on
    legend('Proposed','Location','SE')
end
end
