function val = pitch_of(arr)
    if numel(arr) > 1,
        val = abs(arr(2)-arr(1));
    else
        val = 1;
    end
end
