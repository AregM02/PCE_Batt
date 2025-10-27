% Initializes all necessary data:
% (1) setup a global variable fpath pointing to the desired file location (supported types: .mat, .parquet)
% (2) run the command 'initdata' in the shell
% (!) expected to work flawlessly with metadata from Ahjo/Digatron, adjustments may be needed otherwise;
% (!!) if this script fails, your custom script needs to initialize the following globa variables:
%      - cc/cd: fitted polynomial coefficients for charge/discharge OCV
%      - CNom: nominal capacity
%      - dt: time increment array
%      - soc: state of charge array
%      - t: time array (preferrably in seconds)
%      - vb: voltage array
%      - ib: current array
%      - ecm_Nparams: number of ECM parameters (for example: 5 for 2RC)
%      - ecmfunc: function handle for the ECM @(I,d_t,p,v_init)v2rc(I,d_t,p,v_init)
%      - xrc: initial ECM marameter array (for example: [NaN, NaN, NaN, NaN, NaN] for 2RC)
%      - errors: array for storing errors in parameter estimation (for example: [NaN, NaN, NaN, NaN, NaN] for 2RC)
%      - ix: array to hold the current selected index range (initially has to be 1:length(t))
%      - ix_hyst: index range for the hysteresis test
%      - ix_pulses: index range for the pulse test 
%      - ix_qocv: index range for qocv tests


[~,~,extension] = fileparts(fpath);

if strcmp(extension, '.mat')
    load(fpath)
    data = diga.daten;
    field_names = fieldnames(data);
    for i = 1:length(field_names)
        data.(field_names{i}) = data.(field_names{i})';
    end
    data = rmfield(data, {'T2', 'T3'});
    clear field_names  i diga

elseif strcmp(extension, '.parquet')
    data = parquetread(fpath);
    [~, idx] = unique(data.Programmdauer, 'stable');
    data=data(idx, :);
    data = removevars(data, {'T2', 'T3'});
end
clear extension idx

dt = [diff(data.Programmdauer/1000); 0];
CNom = data.CNom(1);
data.SOC = cumsum(data.Strom.*dt./(3600*CNom));
clear dt

% create necessary global variables
ib = data.Strom; %preprocess current for accurate edge detection
ib = movmean(ib, 10);
ib(abs(ib)<0.005) = 0.;

t = data.Programmdauer/1000;
vb = data.Spannung;
dt = [diff(t); 0];

% define and fix SOC
soc = data.SOC;
idx = find(strcmp(data.Prozedur,'jri_KapTest')); %end of Kaptest: 100% soc
if isempty(idx)
    idx = find(strcmp(data.Prozedur,'jri_KapTest_C2')); % try the other name
end
if ~isempty(idx)
    idx = idx(end);
    soc = soc + 1 - soc(idx);
end
clear idx

% Pulse indices
ix_pulses = find(strcmp(data.Prozedur,'jri_VTC6_PulseHSOC'));

% qOCV indices
ix_qocv = find(strcmp(data.Prozedur,'jri_qOCV_C20'));
ix_qocv = ix_qocv(10:end-10);
ix_qocv_dch = find(ib<-0.005);
ix_qocv_dch = intersect(ix_qocv,ix_qocv_dch);
ix_qocv_cha = find(ib>0.005);
ix_qocv_cha = intersect(ix_qocv,ix_qocv_cha);

% Hysteresis indices
ix_hyst = find(strcmp(data.Prozedur,'jri_hysteresis_socb'));
ix_hyst_cha = ix_hyst(ib(ix_hyst)>0.01);
ix_hyst_dch = ix_hyst(ib(ix_hyst)<-0.01);

% OCV Polynomial
if ~isempty(ix_qocv_cha)
    cd = polyfit(soc(ix_qocv_dch), vb(ix_qocv_dch), 10);
    cc = polyfit(soc(ix_qocv_cha), vb(ix_qocv_cha), 10);
    % c = (cc+cd)/2; %mean ocv
    M = abs(polyval(cc, soc)-polyval(cd, soc))/2; % max hysteresis
end

% ECM function handle for better configurability
ecmfunc = @(I, d_t, p, v_init)v2rc(I, d_t, p, v_init); 
ecm_Nparams = 5;

xrc = NaN(1,ecm_Nparams); % initial ECM parameters
errors = NaN(1,ecm_Nparams); % estimarion errors

%clear data;
% hold on;
% plot(t, vb)
% plot(t, soc)
% hold off;
ix = (1:length(vb))'; %initialize the index
