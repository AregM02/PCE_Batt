% Calculates and plots the simulated voltage over the entire dataset
if exist('c','var') == 1 
    ocv = polyval(c, soc); % if OCV has been fitted (poly coefficients stored in 'c')
else
    reconstocv % otherwise recontruct ocv directly on the voltage profile
    ocv = vb_temp;
    clear vb_temp
end
% v = ocv + vrc_2d(ib, dt, xrc) + hysteresis_2d(ib, dt, M, soc, xh);
v = ocv + ecmfunc(ib, dt, xrc, 0); % simulated voltage

hold on
plot(t(ix), v(ix))
hold off

clear ix_cha;
clear ix_dch;
clear v; % remove this to keep the simulated voltage in memory
clear ocv;

% OUTPUT: none
