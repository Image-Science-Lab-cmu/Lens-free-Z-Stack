function out_val = trianglefilt_and_interp2(xs, ys, vals, xq, yq, x_pitch, y_pitch, xq_pitch, yq_pitch)
    % make triangle filter
    
    vals_filtered = vals;
    if xq_pitch > x_pitch,    
        % from floor(-xq_pitch) to ceil(xq_pitch)
        xx = x_pitch .* (floor(-xq_pitch/x_pitch) : ceil(xq_pitch/x_pitch));
        x_triangle = max(0, 1 - abs(xx)/xq_pitch);
        x_triangle = x_triangle / sum(x_triangle, 'all');
        
        % filter
        vals_filtered = conv2(vals_filtered, reshape(x_triangle, 1, []), 'same');
    end
    
    if yq_pitch > y_pitch,
        % from floor(-yq_pitch) to ceil(yq_pitch)
        yy = y_pitch .* (floor(-yq_pitch/y_pitch) : ceil(yq_pitch/y_pitch));
        y_triangle = max(0, 1 - abs(yy)/yq_pitch);
        y_triangle = y_triangle / sum(y_triangle, 'all');

        % filter
        vals_filtered = conv2(vals_filtered, reshape(y_triangle, [], 1), 'same');
    end

    % interp2
    out_val = interp2(xs, ys, vals_filtered, xq, yq, 'linear', 0);     
end