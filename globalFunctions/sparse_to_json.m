function s = sparse_to_json(A)
%SPARSE_TO_JSON Convert sparse matrix to JSON-friendly struct.

    [rows,cols,vals] = find(A);

    % Convert to zero-based indexing to match Eigen
    rows = rows - 1;
    cols = cols - 1;

    s = struct();

    s.rows = size(A,1);
    s.cols = size(A,2);

    s.trip_rows = rows;
    s.trip_cols = cols;
    s.trip_vals = vals;

end