function A = json_to_sparse(s)
%JSON_TO_SPARSE Convert JSON struct into sparse matrix.

    rows = s.trip_rows + 1;   % convert back to MATLAB indexing
    cols = s.trip_cols + 1;
    vals = s.trip_vals;

    A = sparse(rows,cols,vals,s.rows,s.cols);

end