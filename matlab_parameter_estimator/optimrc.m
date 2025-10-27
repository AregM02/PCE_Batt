% Fits ECM parameters at the given index range 'ix'
% Identifies fall and rise edges in current 'ib'
rise_edges = find(ib(2:end) & ~ib(1:end-1));
rise_edges = intersect(rise_edges, ix);
ib = flip(ib);
fall_edges = find(ib(2:end) & ~ib(1:end-1));
fall_edges = length(ib) - fall_edges;
fall_edges = intersect(fall_edges, ix);
ib = flip(ib);
sort(rise_edges);
sort(fall_edges);

% values are stored as a matrix for all pulses, then averaged
xrc_arr = zeros(numel(rise_edges), 5);
errors_arr = zeros(numel(rise_edges), 5);
errors = zeros(5);

% hold off;
for k=1:length(rise_edges)

    % find voltage jumps
    dv = abs(diff(vb(rise_edges(k):fall_edges(k))));
    ix_jumps = rise_edges(k):fall_edges(k);
    ix_jumps = ix_jumps(dv>0.015);
    
	% check for false identifications
    for j=numel(ix_jumps):-1:2
        if abs(ix_jumps(j-1)-ix_jumps(j))<6 % if too close
            ix_jumps(j-1) = [];
        end
    end
    clear j;
    if isscalar(ix_jumps) || isempty(ix_jumps)
        continue
    end

    % check for CV pulses ans skip over them
    if std(vb(ix_jumps(1)+10:ix_jumps(2)))<0.002
        continue
    end

    % -----Determine R0-------
    dv_r0 = dv(ix_jumps(1)-rise_edges(k)+1);
    xrc(1)=dv_r0/abs(ib(round((fall_edges(k)+rise_edges(k))/2))); % R0
    
    %  ------Fit remaining parameters using the extracted R0------
    width = fall_edges(k) - rise_edges(k); % width of the current pulse
    ixp = rise_edges(k):(rise_edges(k)+10*width);
    % ixp = (fall_edges(i)+10):(fall_edges(i)+4*width);
    ixp = intersect(ixp, ix); % prevent out-of-bounds error

    
    % GENERATE AND REMOVE OCV
    reconstocv; %returns: vb_temp
    vb_temp = vb-vb_temp;
    vb_temp2 = vb_temp; 

    % REMOVE R0*I
    ib_cleaned = zeros(size(ib));
    ib_cleaned(ix_jumps(1)+1:ix_jumps(2)) = ib(round((ix_jumps(1)+ix_jumps(2))/2));
    vb_temp = vb_temp - xrc(1)*ib_cleaned;
    
    % EXPONENTIAL FITTING
    % ====(1) first fit the decay to find tau=RC===========================
        n=3; % index to skip to (safeguard for edge anomalies)
    t_fit = t(ix_jumps(2)+n: ixp(end));
    t_fit = t_fit - t_fit(1);
    v_fit = vb_temp(ix_jumps(2)+n: ixp(end));
    f = fit(t_fit, v_fit, 'exp2');
    tau1 = -1/f.b;
    tau2 = -1/f.d;
    % %----Test Plot----
    % disp(f)
    % plot(t_fit, v_fit)
    % hold on;
    % plot(t_fit, f.a*exp(f.b*t_fit) + f.c*exp(f.d*t_fit))
    % hold off
    % return
    % %-----------------
    % % ---TEST---
    % disp(f)
    % clf; plot(t(ix), vb_temp(ix));
    % hold on; plot(t(ix_jumps(2)+n: ixp(end)), v_fit);
    % plot(t(ix_jumps(2)+n: ixp(end)), f.a*exp(f.b*t_fit) + f.c*exp(f.d*t_fit))
    % return
    % % ----------

    % confidence intervals
    cb = confint(f);
    cb = cb(1,:);
    delta_tau1_inv = abs(cb(2)-f.b);
    delta_tau2_inv = abs(cb(4)-f.d);
    %======================================================================

    % ==== (2) when current is nonzero: R = V/I*(1-exp(1-exp(-t/tau)))=====
    n = 0;
    t_fit = t(ix_jumps(1)+n: ix_jumps(2)-n);
    t_fit = t_fit-t_fit(1);
    v_fit = vb_temp(ix_jumps(1)+n: ix_jumps(2)-n);
    I = ib_cleaned(round((ix_jumps(1) + ix_jumps(2))/2));
    if isempty(v_fit)
        continue
    end
    % R1 = v_fit./(ib_fit.*(1-exp(-t_fit/tau1)));
    ft = fittype('I*R1*(1-exp(-t/tau1)) + I*R2*(1-exp(-t/tau2))', ...
             'independent', 't', 'problem', {'tau1', 'tau2', 'I'});
    fo = fitoptions('Method', 'NonlinearLeastSquares', ...
                'Lower', [0., 0.], 'StartPoint', [1., 1.]);
    f = fit(t_fit, v_fit-v_fit(1), ft, fo, 'problem', {tau1, tau2, I});
    % errors
    cb = confint(f);
    cb = cb(1,:);
    delta_R1 = abs(cb(1)-f.R1);
    delta_R2 = abs(cb(2)-f.R2);
    % calculate the bounds for the relevant variable 1/C with delta method
    delta_C1_inv = sqrt((delta_R1/tau1)^2+(delta_tau1_inv*f.R1)^2);
    delta_C2_inv = sqrt((delta_R2/tau2)^2+(delta_tau2_inv*f.R2)^2);
    % %---TEST---
    % disp(f)
    % clf; plot(t(ix), vb_temp(ix));
    % hold on; plot(t(ix_jumps(1)+n: ix_jumps(2)-n), v_fit);
    % plot(t(ix_jumps(1)+n: ix_jumps(2)-n), v_fit(1)+f.R1*I*(1-exp(-t_fit/tau1)) + f.R2*I*(1-exp(-t_fit/tau2)))
    % hold off
    % return
    % %----------
    %======================================================================

    % ==== R0 CORRECTION (r0+=(v-v_model)/I) ==============================
    % ix_temp = ix_jumps(1)+20:ix_jumps(2)-20;
    % % hold off; plot(vb_temp2(ix));
    % % hold on; plot(ecmfunc(ib(ix_temp), dt(ix_temp), xrc, 0.));
    % % return
    % f_obj = @(x)(sum((vb_temp2(ix_temp) - ecmfunc(ib(ix_temp), dt(ix_temp), [x, xrc(2:end)], 0.) ).^2));
    % x0 = xrc(1);
    % [x0, ~] = fminsearch(f_obj,x0);
    % xrc(1) = x0(1);
    % delta_R0 = 0;

        n=5;
    t_fit = t(ix_jumps(1)+n: ix_jumps(2)-1) - t(ix_jumps(1)+n);
    v_fit = vb(ix_jumps(1)-10) + xrc(1)*I + f.R1*I.*(1-exp(-t_fit/tau1)) + f.R2*I.*(1-exp(-t_fit/tau2));
    dv_r0 = vb(ix_jumps(1)+n: ix_jumps(2)-1) - v_fit;
    dv_r0 = min(abs(dv_r0));
    xrc(1) = xrc(1) + dv_r0/abs(I);
    delta_R0 = dv_r0/abs(I);
    % =====================================================================

    % save results
    xrc = [xrc(1), f.R1, tau1, f.R2, tau2];
    errors = [delta_R0, delta_C1_inv, delta_tau1_inv,...
                        delta_C2_inv, delta_tau2_inv];
    
    % sort xrc
    if xrc(5)<xrc(3)
        temp = xrc;
        temp1 = errors;
        xrc([2 3]) = temp([4 5]);
        xrc([4 5]) = temp([2 3]);

        errors([2 3]) = temp1([4 5]);
        errors([4 5]) = temp1([2 3]);

        clear temp;
    end

    xrc_arr(k,:) = xrc;
    errors_arr(k,:) = errors;


end
clear k dv_r0 dv fall_edges rise_edges width ixp ix_jumps...
    vb_temp vb_temp2 ix_temp f_obj x0 x ib_cleaned t_fit v_fit f I cb ...
    delta_R0 delta_R1 delta_tau1_inv delta_R2 delta_tau2_inv fo ft tau1...
    tau2 delta_C1_inv delta_C2_inv;

if size(xrc)>1
    xrc = mean(xrc_arr);
end

if size(xrc)>1
    errors = mean(errors_arr);
end

% calcv; % show result

% OUTPUT: updated 'xrc' and 'errors' arrays