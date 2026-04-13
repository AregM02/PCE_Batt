% Call this and select a region in the plot with two mouse clicks
% Sets the variable 'ix' to the indices of the selected region
% (!) A plot must be open before running the command

[points, ~] = ginput(2);

% Save the previous selection in case of wanting to go back (command 'b')
if ~exist('ix_prev', 'var')
    ix_prev = {ix};
else
    if length(ix_prev) == 5
        ix_prev = ix_prev(2:end);
    end
    ix_prev{end+1} = ix;
end

ix = find(t>points(1) & t<points(2));
hold off
plot(t(ix), vb(ix))
clear points

% OUTPUT: 'ix_prev' previously selected range
%         'ix' current selected range