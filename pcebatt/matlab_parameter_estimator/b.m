% Goes back to the previous slection in the plot. 
% If no selection has been made or no plot is open, the full dataset will be plotted.
if exist('ix_prev','var') == 1
    if ~isempty(ix_prev)
        hold off
        plot(t(ix_prev{end}), vb(ix_prev{end}))
        ix = ix_prev{end};
        ix_prev = ix_prev(1:end-1);
    else
        hold off
        plot(t, vb)
    end
else
    disp('No previous selection found. Plotting full dataset...')
    plot(t, vb);
end

% OUTPUT: none