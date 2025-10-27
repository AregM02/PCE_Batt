% Finds all pulses in the selected region (within indices 'ix')
% Does not consider major discharge sections (identified through change in SOC)

rising_edges = find(ib(2:end) & ~ib(1:end-1)) + 1;
ixp = rising_edges(rising_edges>=ix(1) & rising_edges<ix(end));

% Set length for the selected pulse region (in # of indices)
if length(ixp)>1
    dist = round((ixp(2)-ixp(1))/9); % 1/27 of the distance between two pulses (change if needed)
else
    dist = ix(end-10) - ixp; % until the end
end

% Check for very large/very small drops in SOC <=> not considered pulses
% Thresholds are 7% and 0.1% respectively
for k = length(ixp):-1:1
    if isscalar(ixp)
        if abs(soc(ixp) - soc(ixp+dist))>0.07 || abs(soc(ixp) - soc(ixp+dist)) < 0.001
            ixp=[];
        end
    elseif abs(soc(ixp(k)) - soc(ixp(k)+2*dist))>0.07 || abs(soc(ixp(k)) - soc(ixp(k)+2*dist))<0.001
        ixp(k)=[];
    end
end
clear k

ixp = ixp + (0:dist);
ixp = ixp(:);
clear rising_edges;
clear dist
ixp = unique(sort(ixp)); % indices of selected regions in matrix form 


%-----------
% Plot
hold on
diff_indices = diff(ixp);
region_starts = ixp([1; find(diff_indices > 1) + 1]);
region_ends = ixp([find(diff_indices > 1); end]);

for j = 1:length(region_starts)
    xregion([t(region_starts(j)), t(region_ends(j))]);
end

hold off

clear j
clear region_ends
clear region_starts
clear diff_indices
%-----------

% OUTPUT: 'ixp' indices of all pulses