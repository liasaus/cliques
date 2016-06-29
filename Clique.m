classdef Clique
    % CLIQUE creates cliques and performs statistical operations on them
    
    properties
        
        N_NODES
        N_CLIQUES
        N_EDGES
        
        % Number of nodes per clique
        CLIQUE_SIZE
        
        % [N_CLIQUES x CLIQUE_SIZE] matrix storing ids of each clique member
        cliques
        
        SIGMOID_SCALE = 24
    end
    
    methods
        
        function C = Clique(n_nodes, n_cliques, clique_size)

            C.N_NODES = n_nodes;
            C.N_CLIQUES = n_cliques;
            C.CLIQUE_SIZE = clique_size;
            C.N_EDGES = round((C.CLIQUE_SIZE^2 / 2) * C.N_CLIQUES);
            
            C.cliques = C.generateCliques();
        end
        
        function cliques = generateCliques(C)
            % Build Cliques matrix (size = [N_CLIQUES x CLIQUE_SIZE] by 
            % randomly drawing integers in the range [1, N_NODES]
            cliques = randi(C.N_NODES,C.N_CLIQUES, C.CLIQUE_SIZE, 'uint32');
            cliques = sort(cliques,2);
        end
        
        function [mean_overlap, overlap_matrix] = findOverlap(C)
        % FINDOVERLAP measures the number of overlapping nodes in each pair
        % of cliques. OVERLAP_MATRIX is a matrix containing each of these 
        % ~N_CLIQUES^2 overlaps. MEAN_OVERLAP is the mean of this matrix.
            overlap_matrix = zeros(C.N_CLIQUES);
            to_compute = tril(ones(C.N_CLIQUES),-1);
            
            for i=1:C.N_CLIQUES
                for j=1:C.N_CLIQUES
                    if to_compute(i,j)
                        % warning: ismembc is an undocumented matlab function!
                        intersection = ismembc(C.cliques(i,:), C.cliques(j,:)); 
                        overlap_matrix(i,j) = sum(intersection);
                    elseif i==j
                        overlap_matrix(i,j) = NaN;
                    end
                end
            end
            % copy lower triangle into upper triangle
            overlap_matrix = overlap_matrix + overlap_matrix';
            mean_overlap = nansum(overlap_matrix(:)) / (C.N_CLIQUES * (C.N_CLIQUES - 1));
        end
        
        function cliques_stimulation = stimulate(C, active_nodes)
        % STIMULATE returns the fraction [0-1] of each clique that is activated
        % when the nodes contained in the vector ACTIVE_NODES in an
        % N_CLIQUES-long vector CLIQUES_STIMULATION
            active_nodes = sort(uint32(active_nodes));
            cliques_stimulation_raw = false(size(C.cliques));
            for i=1:C.N_CLIQUES
                cliques_stimulation_raw(i,:) = ismembc(C.cliques(i,:),active_nodes);
            end
            cliques_stimulation = mean(cliques_stimulation_raw, 2);
        end
        
        function [stimulation, active_nodes] = stimulateRandomNodes(C, p_stimulation)
        % STIMULATERANDOMNODES randomly stimulates with probability
        % P_STIMULATION, and then returns the vector of fractions of clique
        % activations.
            n_stimulated = round(p_stimulation * C.N_NODES);
            active_nodes = randperm(C.N_NODES,n_stimulated); % set of stimulated neurons
            stimulation = C.stimulate(active_nodes);
        end
        
        function firing = stimulation2firing(C, stimulation, threshold)
        % STIMULATION2FIRING converts a vector that denotes the fraction of
        % each clique that is active, STIMULATION, into a binary vector
        % FIRING that indicates whether each clique is fully activated. By
        % default the probability of firing is calculated using a sigmoid
        % transfer function fit to match the one in Miller & Zucker 1999.
        % However, a manually provided THRESHOLD can also be used.
            if nargin > 2
                firing = stimulation >= threshold;
            else
                p_firing = (1 + exp(-C.SIGMOID_SCALE * (stimulation - 0.5))).^-1;
                firing = rand(length(stimulation),1) <= p_firing;
            end
        end
        
        function [clique_error, n_firing_initial, n_firing_final] = predictCliqueError(C, p_input)
                        
            [stimulation, active_nodes] = C.stimulateRandomNodes(p_input);
            cliques_firing = C.stimulation2firing(stimulation);
            
            n_firing_initial = sum(cliques_firing);
            
            cliques_firing_matrix = C.cliques(cliques_firing);
            cliques_firing_nodes = cliques_firing_matrix(:);
            
            % Do the same routine a second time to figure out how many
            % false activations there are
            new_active_nodes = unique([active_nodes'; cliques_firing_nodes]);
            new_stimulation = C.stimulate(new_active_nodes);
            new_cliques_firing = C.stimulation2firing(new_stimulation);
            n_firing_final = sum(new_cliques_firing);
            
            clique_error = (n_firing_final - n_firing_initial) / n_firing_initial;
            
%             fprintf('n_firing_initial = %d\n',n_firing_initial);
%             fprintf('n_firing_final = %d\n',n_firing_final);
%             fprintf('clique_error = %.2f\n',clique_error);
        end
        
        function edges = toEdgeList(C, filename)
        % TOEDGELIST ouptuts the graph in edge list format, either in an
        % [N_EDGES x 2] matrix EDGES or in a .csv file if FILENAME is
        % provided.
        
            edges = zeros(C.N_EDGES * 2,2,'uint16');
            edge_i = 1;

            %% Generate edge list
            for i=1:C.N_CLIQUES
                for j=1:C.CLIQUE_SIZE
                    for k=1:C.CLIQUE_SIZE
                        edges(edge_i,:) = [C.cliques(i,j), C.cliques(i,k)];
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
            
            if nargin == 2
                csvwrite(filename, edges);
            end
        end
        
        function adj = toAdjacencyMatrix(C, filename)
        % TOADJECENCYMATRIX returns ADJ, an adjacency matrix representation
        % of the graph
            adj = zeros(C.N_NODES);

            for i=1:C.N_CLIQUES
                for j=1:C.CLIQUE_SIZE
                    for k=1:C.CLIQUE_SIZE
                        adj(C.cliques(i,j),C.cliques(i,k)) = 1;
                        adj(C.cliques(i,k),C.cliques(i,j)) = 1;
                    end
                end
            end

            % remove diagonal
            for i=1:C.N_NODES
                adj(i,i) = 0;
            end
            
            % outputs to CSV. This file cannot be opened by Gephi!
            % Please use Clique.toEdgeList() if you want to export a graph
            % to Gephi.
            if nargin == 2
                csvwrite(filename, adj);
            end
        end
    end
    
end

