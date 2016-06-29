% GENERATE_GRAPH creates an adjacency matrix representation of a clique
% graph.

N_NEURONS = 10000;
N_CLIQUES = 1000;
CLIQUE_SIZE = 33;
FORMAT = 'edge_list'; % either 'edge_list' or 'adjacency'
FILENAME = 'cliques.csv';

n_synapses = CLIQUE_SIZE^2 * N_CLIQUES;

% Build Cliques
cliques = zeros(N_CLIQUES, CLIQUE_SIZE);
for i=1:N_CLIQUES
    cliques(i,:) = randperm(N_NEURONS,CLIQUE_SIZE);
end

if strcmp(FORMAT,'edge_list')==1
    edges = zeros(n_synapses,2);
    edge_i = 1;
    
    %% Generate edge list
    for i=1:N_CLIQUES
        for j=1:CLIQUE_SIZE
            for k=1:CLIQUE_SIZE
                edges(edge_i,:) = [cliques(i,j), cliques(i,k)];
                edge_i = edge_i + 1;
            end
        end
    end
    
    %% Clean up edge list
    
    % Remove self-weights, e.g. (1,1)
    edges = edges(edges(:,1)~=edges(:,2),:);
    
    % Remove duplicate edges (we must sort first, since [2,3] and [3,2] are
    % the same edge but would not be counted as duplicates
    edges = sort(edges,2);
    edges = unique(edges,'rows');
    
    %% Write to File
    csvwrite(FILENAME, edges);
else
    adj = zeros(N_NEURONS);

    for i=1:N_CLIQUES
        for j=1:CLIQUE_SIZE
            for k=1:CLIQUE_SIZE
                adj(cliques(i,j),cliques(i,k)) = 1;
                adj(cliques(i,k),cliques(i,j)) = 1;
            end
        end
    end

    % remove diagonal
    for i=1:N_NEURONS
        adj(i,i) = 0;
    end

    adj_gephi = zeros(N_NEURONS+1);
    adj_gephi(2:N_NEURONS+1,2:N_NEURONS+1) = adj;
    adj_gephi(1,2:N_NEURONS+1) = 1:N_NEURONS;
    adj_gephi(2:N_NEURONS+1,1) = 1:N_NEURONS;
end
