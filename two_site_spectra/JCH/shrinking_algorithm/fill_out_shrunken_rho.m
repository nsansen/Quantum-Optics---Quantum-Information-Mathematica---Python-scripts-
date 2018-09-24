function filled_rho = fill_out_shrunken_rho(full_dim, min_evector, connected_elements_store, non_repeated_elements)

filled_rho = sparse(full_dim, full_dim);

% Fill in the subset of values that were explicitly solved for:
filled_rho(non_repeated_elements) = min_evector;
num_connected_paths = length(connected_elements_store);

for loop = 1:num_connected_paths
    
    if ne(length(connected_elements_store{loop}),1) % If this path is not just an element connected to itself
       
        % Add the columns of L corresponding to the redundant elements:
        elements_to_fill_this_connection = connected_elements_store{loop}(2:end);
        filled_rho(elements_to_fill_this_connection) = filled_rho(connected_elements_store{loop}(1));
        
    else
        
        continue
        
    end
    
end

filled_rho = filled_rho/trace(filled_rho);

end