function matrix_big = tensor_matrix(matrix, M, index)
% Function to return the tensor product: I * I * ... * matrix * I * ... *
% I, where matrix is in the index'th position, and '*' denotes a tensor
% multiplication, and there are M matrices all together in the product

id = speye(size(matrix)); % The identity matrix to be tensored successively
matrix_big = 1; % Initialise the output matrix (NB kron(1, matrix) = matrix

for k = 1:M
    
    if (k == index)
        next_matrix = matrix;
    else
        next_matrix = id;
    end
    
    matrix_big = kron(matrix_big, next_matrix);
    
end

end
