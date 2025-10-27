% Records parameter values 'xrc' to a table named 'results'
% Requires: xrc (ECM parameter array), errors (errors in their estimation)

soc_temp = 100*mean(soc(ix)); % avg soc
temperature = mean(data.Temp(ix));

% Snap value to the nearest grid point (optional)
%grid = 0:5:100;
%soc_temp = interp1(grid, grid, soc_temp, 'nearest', 'extrap');

current_row = table(xrc(1), xrc(2), xrc(3), xrc(4), xrc(5), ...
                    errors(1), errors(2), errors(3), errors(4), errors(5), ...
                    soc_temp, temperature, ...
                    VariableNames={'R0', 'R1', 'tau1', 'R2', 'tau2', ...
                              'dR0', 'dC1_inv', 'dtau1_inv', 'dC2_inv', 'dtau2_inv', ...
                              'SOC', 'T'});

if exist('results','var') == 1
    results = [results; current_row];
else
    results = current_row;
end
clear current_row soc_temp grid temperature;

% OUTPUT: 'results' tabe containing xrc and errors at given temperature and SOC
