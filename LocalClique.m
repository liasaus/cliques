classdef LocalClique < Clique
    % LOCALCLIQUE creates cliques that are wired together based on the 
    % distances between them, creating local structure.
    
    properties
        
    end
    
    methods
        
        function LC = LocalClique(n_nodes, n_cliques, clique_size)
            LC = LC@Clique(n_nodes, n_cliques, clique_size);
        end
        
        function cliques = generateCliques(C)
            % Build CLIQUES matrix (size = [N_CLIQUES x CLIQUE_SIZE])
            
            % @todo wire based on distance rule
            cliques = randi(C.N_NODES,C.N_CLIQUES, C.CLIQUE_SIZE, 'uint32');
            cliques = sort(cliques,2);
        end
       
    end
    
end

