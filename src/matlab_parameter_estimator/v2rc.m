% 2-RC ECM
function v_total = v2rc(strom, dt, params, v_0)
    % Calculates output of a 2RC ECM
    % strom: current
    % dt: time increment array
    % params: currently expected input [r0, r1, tau1, r2, tau2]
    % v_0: initial voltage

    % Parse parameters
    R0 = params(1);
    R1 = params(2);
    C1 = params(3)/R1;
    R2 = params(4);
    C2 = params(5)/R2;

    % Number of time steps
    N = size(strom);
    
    % Precompute a and b arrays for the first RC pair
    a1 = exp(-dt ./ (R1 * C1));
    b1 = strom .* (1 - exp(-dt ./ (R1 * C1))) .* R1;

    % Precompute a and b arrays for the second RC pair
    a2 = exp(-dt ./ (R2 * C2));
    b2 = strom .* (1 - exp(-dt ./ (R2 * C2))) .* R2;

    % Initialize recursive voltage vectors for both RC pairs
    v_temp1 = zeros(size(strom));
    v_temp2 = zeros(size(strom));

    % Initial conditions
    v_temp1(1) = b1(1);
    v_temp2(1) = b2(1);

    % Compute voltages for the RC pairs recursively
    for j = 2:numel(strom)
        v_temp1(j) = a1(j) * v_temp1(j-1) + b1(j);
        v_temp2(j) = a2(j) * v_temp2(j-1) + b2(j);
    end
   
    % Total voltage = ohmic + RC contributions
    v_total = R0 * strom + v_temp1 + v_temp2;

    % Add initial voltage, if provided 
    v_total = v_total + v_0;
    
    % If params are invalid (usually after initialization), return zeros
    if any(isnan(params(1:5)))
        v_total = zeros(N);
    end
end