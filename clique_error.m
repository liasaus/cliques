CLIQUE_SIZE = 30;
CLIQUE_DENSITIES = [0.1, 0.5, 1, 2, 3];
NODE_COUNTS = [10000, 50000, 100000, 200000, 300000 ];
N_DRAWS = 20;
P_INPUT = 0.2;

error_matrix = zeros(length(CLIQUE_DENSITIES),length(NODE_COUNTS));
initial_matrix = zeros(length(CLIQUE_DENSITIES),length(NODE_COUNTS));
final_matrix = zeros(length(CLIQUE_DENSITIES),length(NODE_COUNTS));
for i=1:length(CLIQUE_DENSITIES)
    error_js = zeros(length(NODE_COUNTS),1);
    initial_js = zeros(length(NODE_COUNTS),1);
    final_js = zeros(length(NODE_COUNTS),1);
    for j=1:length(NODE_COUNTS)
        n_cliques = round(NODE_COUNTS(j)*CLIQUE_DENSITIES(i));
        fprintf('%d nodes, %d cliques\n',NODE_COUNTS(j),n_cliques);
        errors = zeros(N_DRAWS,1);
        initials = zeros(N_DRAWS, 1);
        finals = zeros(N_DRAWS, 1);
        parfor k=1:N_DRAWS
            C = Clique(NODE_COUNTS(j), n_cliques, CLIQUE_SIZE);
            [errors(k), initials(k), finals(k)] = C.predictCliqueError(P_INPUT);
        end
        error_js(j) = nanmean(errors);
        initial_js(j) = nanmean(initials) / n_cliques;
        final_js(j) = nanmean(finals) / n_cliques;
    end
    error_matrix(i,:) = error_js;
    initial_matrix(i,:) = initial_js;
    final_matrix(i,:) = final_js;
end

figure();
imagesc(error_matrix);
colorbar();

ylabels = cellstr(num2str(CLIQUE_DENSITIES', '%.1f'));
set(gca,'YTick',1:length(CLIQUE_DENSITIES));
set(gca,'YTickLabel',ylabels);
set(gca,'YDir','normal');

xlabels = cellstr(num2str(NODE_COUNTS', '%d'));
set(gca,'XTick',1:length(NODE_COUNTS));
set(gca,'XTickLabel',xlabels);
xlabel('nodes');
ylabel('clique density');
title('Error');


% initial and final matrices should use same color scale
color_scale = [min([initial_matrix(:);final_matrix(:)]), max([initial_matrix(:);final_matrix(:)])];

figure();
imagesc(initial_matrix);
colorbar();

ylabels = cellstr(num2str(CLIQUE_DENSITIES', '%.1f'));
set(gca,'YTick',1:length(CLIQUE_DENSITIES));
set(gca,'YTickLabel',ylabels);
set(gca,'YDir','normal');

xlabels = cellstr(num2str(NODE_COUNTS', '%d'));
set(gca,'XTick',1:length(NODE_COUNTS));
set(gca,'XTickLabel',xlabels);
xlabel('nodes');
ylabel('clique density');
title('Fraction of Active Cliques (Initial)');
caxis(color_scale);

figure();
imagesc(final_matrix);
colorbar();

ylabels = cellstr(num2str(CLIQUE_DENSITIES', '%.1f'));
set(gca,'YTick',1:length(CLIQUE_DENSITIES));
set(gca,'YTickLabel',ylabels);
set(gca,'YDir','normal');

xlabels = cellstr(num2str(NODE_COUNTS', '%d'));
set(gca,'XTick',1:length(NODE_COUNTS));
set(gca,'XTickLabel',xlabels);
xlabel('nodes');
ylabel('clique density');
title('Fraction of Active Cliques (Final)');
caxis(color_scale);