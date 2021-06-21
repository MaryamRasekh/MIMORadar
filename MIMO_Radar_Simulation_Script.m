close all;
clear;
clc;
single_run = true;
monte_carlo_run = false;

SimParameters = struct(); % if empty will use some default parameters
SimParameters.N_tx = 128; % number of elements in TX phased array
% SimParameters.N_rx = 1; % number of elements in digital receiver
SimParameters.N_beacons = 15; % number of compressive beacons (subframes)
SimParameters.N_chirp = 32; % number of chirps in each subframe
SimParameters.N_symb = 256;  % number of samples in a single chirp
SimParameters.perSymb_SNR_dB = -5; % mean per sample SNR when one transmitter is active
SimParameters.T_gap = 512; % duration between consecutive chirps (as multiples of symbol period)
SimParameters.DR = 10; % dynamic range of target signal amplitudes

%% single run for visual representation
if single_run
    PlotResultsFlag = true; % plot the pretty pictures
    N_target = 6;
    SimulateMIMORadarFrame(N_target,SimParameters,PlotResultsFlag);
end

%% run monte carlo to get error cdf
if monte_carlo_run
    N_tot_targets = 20; % total number of targets to be simulated
    min_targets = 5; % minimum number of targets in each realization 
    max_targets = 8; % maximum number of targets in each realization
    PlotResultsFlag = false; % don't plot the pretty pictures
    error_mat_master = zeros(N_tot_targets,3); % contains absolute errors in range, doppler, spatial frequency 
                                               % (in units of corresponding FFT grid size)
    number_of_simulated_targets = 0;
    while number_of_simulated_targets<N_tot_targets
        if N_tot_targets - number_of_simulated_targets <= max_targets
            N_target = N_tot_targets - number_of_simulated_targets;
        else
            N_target = min_targets + round(rand()*(max_targets-min_targets));
        end
        error_mat = SimulateMIMORadarFrame(N_target,SimParameters,PlotResultsFlag);
        error_mat_master(number_of_simulated_targets+1:number_of_simulated_targets+N_target,:) = error_mat;
        number_of_simulated_targets = number_of_simulated_targets + N_target;
        clc;  disp([num2str(number_of_simulated_targets),' out of ',num2str(N_tot_targets),' done.'])
    end
    figure()
    subplot(3,1,1)
    cdfplot(error_mat_master(:,1)); title('range error (normalized to FFT grid size)')
    subplot(3,1,2)
    cdfplot(error_mat_master(:,2)); title('doppler error (normalized to FFT grid size)')
    subplot(3,1,3)
    cdfplot(error_mat_master(:,3)); title('direction error (normalized to FFT grid size)')
end


% --------------------------------------------------------------------------------------------------------------
%% main function that runs sim:
% --------------------------------------------------------------------------------------------------------------

function error_mat = SimulateMIMORadarFrame(N_target,SimParameters,PlotResultsFlag)
% Simulates one frame of target acquisition for settings provided in
% SimParameters object. If none or incorrectly provided, default settings
% used.
% N_target: number of targets present in the scene
% PlotResultsFlag: if true, will plot range-doppler heatmaps after each
% target extraction

%% initialize and simulate frame response
%%% Fields required in SimParameters object:
%%% (initialized with default values)
N_tx = 128; % number of elements in TX phased array
%   N_rx = 1; % number of elements in digital receiver
N_beacons = 15; % number of compressive beacons (subframes)
N_chirp = 32; % number of chirps in each subframe
N_symb = 256;  % number of samples in a single chirp
perSymb_SNR_dB = -5; % mean per sample SNR when one transmitter is active
T_gap = 512; % duration between consecutive chirps (as multiples of symbol period)
DR = 10; % dynamic range of target signal amplitudes

%%% override parameters with input values if provided
if isfield(SimParameters, 'N_tx')
    N_tx = SimParameters.N_tx;
end
% if isfield(SimParameters, 'N_rx')
%     N_rx = SimParameters.N_rx;
% end
if isfield(SimParameters, 'N_beacons')
    N_beacons = SimParameters.N_beacons;
end
if isfield(SimParameters, 'N_chirp')
    N_chirp = SimParameters.N_chirp;
end
if isfield(SimParameters, 'N_symb')
    N_symb = SimParameters.N_symb;
end
if isfield(SimParameters, 'perSymb_SNR_dB')
    perSymb_SNR_dB = SimParameters.perSymb_SNR_dB;
end
if isfield(SimParameters, 'T_gap')
    T_gap = SimParameters.T_gap;
end
if isfield(SimParameters, 'DR')
    DR = SimParameters.DR;
end
noise_power = 10^(-perSymb_SNR_dB/10);

%%% fast FFT frequency (proportional to delay)
ffreq_guard_bins = min(2,floor(N_symb*0.01));
min_ffreq = ffreq_guard_bins*2*pi/N_symb;
max_ffreq = 2*pi - ffreq_guard_bins*2*pi/N_symb;

%%% doppler offsets, unit: radians per chirp duration (total time between start of two consecutive chirps)
doppler_guard_bins = min(2,floor(N_chirp*0.01));
min_doppler = -pi + doppler_guard_bins*2*pi/N_chirp;
max_doppler =  pi - doppler_guard_bins*2*pi/N_chirp;

%%% angles, represented by SPATIAL FREQUENCY on the (linear) TX array
omega_guard_bins = min(2,floor(N_tx*0.01));
min_omega = -pi + omega_guard_bins*2*pi/N_tx;
max_omega =  pi - omega_guard_bins*2*pi/N_tx;

%%% randomly draw target delay, doppler, angles (equivalent spatial frequencies)
ffreq_vect = min_ffreq + (max_ffreq-min_ffreq)*rand(1,N_target); % maybe change this to agree with uniform distribution in area (likelihood increases linearly with delay)
doppler_vect = min_doppler + (max_doppler - min_doppler)*rand(1,N_target);
omega_vect = min_omega + (max_omega-min_omega)*rand(1,N_target);

%%% calculate complex amplitude of signal received from each target
% This part can probably be done better based on physical models :D
%
%      crude modeling of variation in radar cross section
%           ------------------------
amp_vect = (0.8+0.4*rand(1,N_target))./ffreq_vect.^2;
%                                    ---------------
%                              amplitude fall-off with range
%                   (power attenuation prop. to fourth power of distance)
%
%%% normalize and ensure dynamic range does not exceed limit (DR)
amp_vect = amp_vect./sqrt(mean(amp_vect.^2));
while max(amp_vect)/min(amp_vect) > DR
    amp_vect = amp_vect + 0.1;
    amp_vect = amp_vect./sqrt(mean(amp_vect.^2));
end
amp_vect = amp_vect.*exp(1i*2*pi*rand(1,N_target)); % uniformly random phases for target amplitudes

%%% generate compressive beacon weights
beacon_pool = 1/sqrt(N_tx)*exp(1i*2*pi*rand(N_tx,N_beacons));
A_CS = beacon_pool.';

%%% compute compressive beacon gain for each target direction
beacon_gains = exp(1i*omega_vect'*(0:N_tx-1))*beacon_pool;

%%% generate raw received signal over one frame
y_rx = zeros(N_beacons,N_chirp,N_symb);
for target = 1:N_target
    f_delay = ffreq_vect(target);
    y_nom_chirp = amp_vect(target)*exp(1i*f_delay*(0:N_symb-1)); % nominal chirp response for this target
    dop_phase = exp(1i*doppler_vect(target)*(0:N_chirp-1)'); % phase change over chirps due to doppler (in one subframe)
    y_nom = dop_phase*y_nom_chirp; % nominal subframe response for this target
    beac_dop_phase = exp(1i*doppler_vect(target)*N_chirp*(0:N_beacons-1)); % phase change over subframes due to doppler (in one frame)
    for beacon = 1:N_beacons
        % apply beacon and doppler modulation on subframes
        y_tot_target = beac_dop_phase(beacon)*beacon_gains(target,beacon)*reshape(y_nom,1,N_chirp,N_symb);
        % add this target's signal to the overall received signal
        y_rx(beacon,:,:) = y_rx(beacon,:,:) + y_tot_target;
    end
end
% rx_sig = y_rx; % store noiseless version
%%% add noise
y_rx = y_rx + sqrt(noise_power/2)*(randn(size(y_rx)) + 1i*randn(size(y_rx)));
% Ok, we have our simulated measurements!
% now, let's detect

%% estimate targets from simulated measurements
oversampling_symb = 8; % range FFT oversampling rate
oversampling_chirp = 2^ceil(log(N_beacons)/log(2)+1); % doppler FFT oversampling rate
%%% evaluate range-doppler domain signal
y_rd = zeros(N_beacons, N_chirp*oversampling_chirp, N_symb*oversampling_symb);
for i_beacon = 1:N_beacons
    y_rd(i_beacon,:,:) = reshape(fft2(reshape(y_rx(i_beacon,:,:),N_chirp,N_symb),N_chirp*oversampling_chirp,N_symb*oversampling_symb),1,N_chirp*oversampling_chirp,N_symb*oversampling_symb);
end
% aggregate power in each range-doppler bin from all beacons
power_bins = reshape(sum(abs(y_rd).^2,1),N_chirp*oversampling_chirp,N_symb*oversampling_symb);
peak_power = max(max(power_bins));
if PlotResultsFlag
    % let's choose some parameters to help with units for range and doppler
    % (radial speed)
    c = 3e8; % speed of light
    fc = 60e9; lambda = c/fc; % carrier frequency, wavelength
    BW = 60e6; ts = 1/BW; % bandwidth, symbol period
    chirp_slope = 2*pi*1e13; % slope of chirp frequency ramp in Hz/s
    ffreq_to_range = 1/ts/chirp_slope*c/2; % scaling that translates delay in symbol periods to range in meters
    doppler_to_speed = lambda/2/pi/(ts*(N_symb+T_gap)); % scaling that translates doppler in radians per symbol period to radial speed in meters per second
    % oversampled range-doppler axes:
    range_axis = (0:N_symb*oversampling_symb-1)*2*pi/N_symb/oversampling_symb*ffreq_to_range;
    speed_axis = angle(exp(1i*(0:N_chirp*oversampling_chirp-1)/oversampling_chirp*2*pi/N_chirp))*doppler_to_speed;
    
    min_val = max(max(power_bins))/20000;
    fig_0 = figure();
    %     plot(delay_vect*oversampling_symb+1,(doppler_vect+2*pi*(doppler_vect<0))*N_chirp*oversampling_chirp/2/pi,'rx','LineWidth',.75);
    %     xlim([0,N_symb*oversampling_symb]); ylim([0,N_chirp*oversampling_chirp])
    plot(ffreq_vect*ffreq_to_range,doppler_vect*doppler_to_speed,'rx','LineWidth',.75);
    hold on;
    xlim([0,2*pi*ffreq_to_range]); ylim([-1,1]*pi*doppler_to_speed)
    xlabel('range (m)'); ylabel('radial speed (m/s)');
    title('true (x) and estimated (o) targets')
    fig_1 = figure();
    h = surf(range_axis,speed_axis,20*log10(power_bins));
    set(h,'edgecolor','none'); view(2);
    xlim([0,2*pi*ffreq_to_range]); ylim([-1,1]*pi*doppler_to_speed)
    xlabel('range (m)'); ylabel('radial speed (m/s)');
    title('aggregate bin power')
end

%%% find the (N_target + extra_bins) strongest components in range-doppler
%%% signal. This part is a crude version of NOMP crude and will be cleaned
%%% up in next version.
extra_bins = floor(log(N_target));
target_bins = zeros(N_target+extra_bins,2);
power_residue = power_bins; % updated after each peak is identified and subtracted
y_residue = y_rx; % updated after each peak is identified and subtracted
y_rd_residue = y_rd; % updated after each peak is identified and subtracted

if PlotResultsFlag
    num_subplot_rows = 2; % floor(sqrt(N_target+extra_bins));
    num_subplot_cols = 2; % ceil(N_target+extra_bins/num_subplot_rows);
    fig_2 = figure();
    subplot(num_subplot_rows,num_subplot_cols,1);
end
for i_target = 1:N_target+extra_bins
    [max_v,ind_row_vect] = max(power_residue);
    [~,ind_col] = max(max_v);
    ind_row = ind_row_vect(ind_col);
    target_bins(i_target,:) = [ind_row,ind_col]; % row and column index of largest peak in range-doppler signal
    w_range = 2*pi*(ind_col-1)/N_symb/oversampling_symb; % equivalent range frequency of largest peak
    w_doppler = angle(exp(1i*2*pi*(ind_row-1)/N_chirp/oversampling_chirp)); % equivalent doppler frequency of largest peak (identified target)
    y0 = reshape(exp(1i*w_doppler*(0:N_chirp-1)')*exp(1i*w_range*(0:N_symb-1)),1,N_chirp,N_symb); % time domain response of identified target
    for i_beacon = 1:N_beacons % extract the identified target's response from signal residue in fast-slow time and range-doppler domain
        y_residue(i_beacon,:,:) = y_residue(i_beacon,:,:) - ...
            y_rd_residue(i_beacon,ind_row,ind_col)/N_symb/N_chirp*y0;
        y_rd_residue(i_beacon,:,:) = reshape(fft2(reshape(y_residue(i_beacon,:,:),N_chirp,N_symb),N_chirp*oversampling_chirp,N_symb*oversampling_symb),1,N_chirp*oversampling_chirp,N_symb*oversampling_symb);
    end
    power_residue = reshape(sum(abs(y_rd_residue).^2,1),N_chirp*oversampling_chirp,N_symb*oversampling_symb); % update aggregate power residue
    if PlotResultsFlag
        if i_target <= num_subplot_rows*num_subplot_cols
            figure(fig_2); subplot(num_subplot_rows,num_subplot_cols,i_target);
            h = surf(range_axis,speed_axis,10*log10(power_residue+min_val));
            set(h,'edgecolor','none'); view(2);
            % xlim([0,N_symb*oversampling_symb]); ylim([0,N_chirp*oversampling_chirp]);
            xlim([0,2*pi*ffreq_to_range]); ylim([-1,1]*pi*doppler_to_speed)
            ylabel(['after ',num2str(i_target),' extraction(s)']);
        end
        figure(fig_0); hold on;
        plot((ind_col-1)*2*pi/N_symb/oversampling_symb*ffreq_to_range,w_doppler*doppler_to_speed,'bo','LineWidth',1,'MarkerSize',5+5*power_bins(ind_row,ind_col)/peak_power)
    end
end

%%% now find the direction of each target (or targets inside each
%%% identified bin, for now just assuming one target in each bin)
omega_est = zeros(N_target,1);
for i_bin = 1:length(target_bins(:,1))
    y_b = y_rd(:,target_bins(i_bin,1),target_bins(i_bin,2));
    y_b = y_b(:);
    dop_i = (target_bins(i_bin,1)-1)*2*pi/N_chirp/oversampling_chirp;
    y_b = y_b.*exp(-1i*dop_i*N_chirp*(0:N_beacons-1)');
    tau = norm(y_b)^2/7;
    [om_list, g_list, r_list] = extractSpectrum_MR(y_b, A_CS, tau, 8, 6);
    omega_est(i_bin) = angle(exp(1i*om_list(1)));
end

%% find absolute error in range, doppler and spatial frequency estimates (in grid size)
% (fine tuned doppler estimation not coded yet)
detected_targets = zeros(N_target+extra_bins,3); % delay, doppler (coarse), spatial frequency
detected_targets(:,1) = (target_bins(:,2)-1)/oversampling_symb; % delay
detected_targets(:,2) = N_chirp/2/pi*angle(exp(1i*2*pi*(target_bins(:,1)-1)/N_chirp/oversampling_chirp));
detected_targets(:,3) = omega_est*N_tx/2/pi;

true_targets = [ffreq_vect(:)*N_symb/2/pi, doppler_vect(:)*N_chirp/2/pi, omega_vect(:)*N_tx/2/pi];
error_mat = zeros(N_target,3);
% for each true target, look for closest match in estimated targets and
% compute error
for i_target = 1:N_target
    error_vects = abs((detected_targets - true_targets(i_target,:)));
    [~,ind_match] = min(sum(error_vects,2));
    error_mat(i_target,:) = error_vects(ind_match,:);
end
end


% --------------------------------------------------------------------------------------------------------------
% --------------------------------------------------------------------------------------------------------------
% --------------------------------------------------------------------------------------------------------------


%% NOMP functions:

function [omegaList, gainList, residueList] = extractSpectrum_MR(y, S,...
			      	   tau, overSamplingRate, numRefine)
% SUMMARY:
% 
%   given measurements: y = S * (mixture of sinusoids) + white noise
%          returns two **ordered** lists: a list of 
%          estimates of frequencies of the sinusoids in the mixture 
%          and list of corresponding gains 
% INPUT:
%    y - measurements
%    S - measurement matrix
%    tau - minimum payoff per sinusoid
%        - should be slightly more than the noise level sigma2
%        minimize: norm(y - sum of sinusoids)^2 + tau * ell0_norm(gain)
%        ** LARGE VALUES OF TAU PROMOTE SPARSITY **% 
%    overSamplingRate (optional) -  # grid points = overSamplingRate * N
%              
%    numRefine (optional) : number of refinement steps after each frequency
%                           
% OUTPUT:
%   omegaList    - frequencies
%   gainList     - gains of estimated frequencies
%   residueList  - trajectory of the energy in the residual

if ~exist('overSamplingRate','var'), overSamplingRate = 4;
elseif isempty(overSamplingRate), overSamplingRate = 4; end

if ~exist('numRefine','var'), numRefine = 6;
elseif isempty(numRefine),    numRefine = 6; end

% Algorithm preprocessing
sampledManifold = preProcessMeasMat(S, overSamplingRate);

if sampledManifold.is_eye
	S = []; % erase identity matrix
end

omegaList = [];
gainList  = [];
y_r = y;

residueList = [ y_r' * y_r ];

while true
    
    % keep detecting new sinusoids until power in residue 
    % becomes small; *** how small *** determined by *** tau ***
    
    % detect gain and frequency of an additional sinusoid
    [omega_new, gain_new, y_r, res_inf_normSq_rot] = ...
        detectNew(y_r, sampledManifold);
    % detecttNew removes the contribution of the newly detected
    % from the old residue  y_r(input) and reports the new residual
    % measurement y_r (output)
    
    % stopping criterion:
    if res_inf_normSq_rot < tau
        break;
    end
    
    % newly detected sinusoid is coarse - so we refine it to 
    % imitate detection on the continuum
    for i = 1:numRefine
        [omega_new, gain_new, y_r] = refineOne(y_r, omega_new, ...
            gain_new, S, sampledManifold.ant_idx, true);
    end
    % refineOne checks whether the refinement step decreases
    % the l-2 norm of the residue y_r
    
    % Add newly detected sinusoid to the ordered lists
    omegaList = [omegaList; omega_new];
    gainList  = [gainList; gain_new];
    
    % refine all frequencies detected so far
    % can be interpreted as a search for better frequency supports
    [omegaList, gainList, y_r] = refineAll(y_r, omegaList,...
        gainList, S, sampledManifold.ant_idx, numRefine);
    % refineAll only uses refineOne to tweak parameters and the energy 
    % in the residual measurements y_r can only decrease as a result

    % Solve least squares for the dictionary set [Ax(omega)] omega in 
    % omegaList
    if false
    [omegaList, gainList, y_r] = solveLeastSquares(y , omegaList, ...
        S, sampledManifold.ant_idx);    
    end
    % ensures that for the support we have settled on our choice of 
    % gains is optimal (in terms giving us the lowest residual energy)
    % though the practical value of the least-squares step is debatable
    % when the frequency support is well-condition (since least-squares
    % is a by-product of refineAll once the frequencies have converged)
    % we need this step for theoretical guarantees on convergence rates
 
    residue_new = y_r'*y_r;
    residueList = [residueList; residue_new];
    
end

% revert to standard notion of sinusoid: 
%           exp(1j*(0:(N-1))'*omega)/sqrt(N)
gainList = gainList .* exp(1j*sampledManifold.ant_idx(1)*omegaList);

end

% --------------------------------------------------------------------

function [omega, gain, y_r, res_inf_normSq_rot] = detectNew(y,...
					 sampledManifold)
% SUMMARY:
% 
% 	detects a new sinusoid on the coarse grid
% 	refines the newly detected frequency *numRefine* times
% 
% INPUT:
% 	y - measurements
% 	sampledManifold - when measurements are compressive,
% 		  we precompute (using IFFT operation) a
% 		  **dictionary** of responses for sinusoids 
%  		  corresponding to frequencies on a coarse grid 
% 		  and store them in the MATLAB **structure** 
% 		  sampledManifold
%
% OUTPUT:
% 	omega - frequency on [0,2*pi) which best explains
% 			the measurements y
% 	gain  - corresponding complex gain
% 	y_r   - after removing the detected sinusoid from the
%         	measurements y, y_r is the ***residual measurement***
%          res_inf_normSq_rot - max energy among DFT directions - needed
%          for stopping criterion

R = length(sampledManifold.coarseOmega);
N = sampledManifold.length;
OSR = round(R/N);

if sampledManifold.is_eye
    gains  = fft(y, R)/sqrt(N);
    
    if sampledManifold.ant_idx(1)~=0
        gains = gains.*exp(-1j*sampledManifold.coarseOmega(:)*...
            sampledManifold.ant_idx(1));
    end
    prob = (abs(gains).^2);
else
    energy = sampledManifold.map_IfftMat_norm_sq.';
    gains = (sampledManifold.map_IfftMat'*y)./energy;
    prob = (abs(gains).^2).*energy;
end

[~,IDX] = max(prob);

omega = sampledManifold.coarseOmega(IDX);
gain = gains(IDX);

% compute the response corresponding to the
% current estimate of the sinusoid
if sampledManifold.is_eye
    x = exp(1j*sampledManifold.ant_idx * omega)...
        /sqrt(sampledManifold.length);
else
    x = sampledManifold.map_IfftMat(:,IDX);
end

% residual measurements after subtracting
% out the newly detected sinusoid
y_r = y - gain * x;

% For stopping criterion
% we check only DFT frequencies - gives us a handle
% the false alarm rate (a measure of how often we would 
% over estimate the size of the support set)

res_inf_normSq_rot = max(prob(1:OSR:end));

end

% --------------------------------------------------------------------

function [omega, gain, y_r] = refineOne(y_r, omega, gain, S,...
			 ant_idx, isOrth)
% SUMMARY:
%   Refines parameters (gain and frequency) of a single sinusoid
%   and updates the residual measurement vector y_r to reflect
%   the refinement
% INPUT:
% 	y_r - residual measurement (all detected sinusoids removed)
%	omega - current estimate of frequency of sinusoid we want to
% 			refine
%	gain - current estimate of gain of sinusoid we want to refine
% 	S - measurement matrix
%           if [], measurements are direct i.e S = eye(N),
% 	    where N is length of sinusoid	
% 	ant_idx - translating indexes to phases in definition of sinusoid
%	isOrth - binary flag - is y_r orthogonal to x(omega) 
%	       - default - false
% OUTPUT:
%       refined versions of omega, gain and y_r
%       (see INPUT for definitions)

if ~exist('isOrth', 'var'), isOrth = false; end

if isempty(S)
    is_eye = true;
else
    is_eye = false;
end

N = length(ant_idx);
x_theta  = exp(1j*ant_idx*omega)/sqrt(N);
dx_theta = 1j * ant_idx .* x_theta;
d2x_theta = -(ant_idx.^2) .* x_theta;

if ~is_eye
    x_theta   = S * x_theta;
    dx_theta  = S * dx_theta;
    d2x_theta = S * d2x_theta;
end

% add the current estimate of the sinusoid to residue
y = y_r + gain*x_theta;


% UPDATE GAIN
% recompute gain and residue to ensure that 
% y_r is orthogonal to x_theta - this property
% is lost when we refine other sinusoids
if ~isOrth
    if is_eye
        energy = 1;
    else
        energy = x_theta'*x_theta;
    end
    gain = (x_theta'*y)/energy;
    y_r = y - gain*x_theta;
end

der1 = -2*real(gain * y_r'*dx_theta);
der2 = -2*real(gain * y_r'*d2x_theta) +...
    2*abs(gain)^2*(dx_theta'*dx_theta);

% UPDATE OMEGA
if der2 > 0
    omega_next = omega - der1/der2;
else
    omega_next = omega - sign(der1)*(1/4)*(2*pi/N)*rand(1);
end

% COMPUTE x_theta for omega_next so that we can compute 
% gains_next and y_r_next
x_theta  = exp(1j*ant_idx*omega_next)/sqrt(N);
if is_eye
    energy = 1;
else
    x_theta = S * x_theta;
    energy = (x_theta'*x_theta);
end

% UPDATE GAIN
gain_next = (x_theta'*y)/energy;

% UPDATE RESIDUE
y_r_next = y - gain_next*x_theta;

% check for decrease in residue -  needed as a result of 
% non-convexity of residue (even when the cost surface 
% is locally convex); This is the same as checking whether 
% |<y, x_theta>|^2/<x_theta,x_theta> improves as we ensured 
% that y - gain*x_theta is perp to x_theta by recomputing gain

if (y_r_next'*y_r_next) <= (y_r'*y_r)
    % commit only if the residue decreases
    omega = omega_next;
    gain = gain_next;
    y_r = y_r_next;
end

end

% --------------------------------------------------------------------

function [omegaList, gainList, y_r] = solveLeastSquares(y , omegaList, ...
    S, ant_idx)
% SUMMARY:
%    Reestimates all gains wholesale
N = length(ant_idx);
if isempty(S) % is this an identity matrix 
    A = exp(1j*ant_idx*omegaList.')/sqrt(N);
else
    A = S * exp(1j*ant_idx*omegaList.')/sqrt(N);
end   

% update gains
gainList = (A'*A)\(A'*y);
% update residues
y_r = y - A*gainList;

% energy in the residual measurement y_r is guaranteed to not increase 
% as a result of this operation. Therefore, we refrain from checking 
% whether the energy in y_r has increased
end

% --------------------------------------------------------------------

function [omegaList, gainList, y_r] = refineAll(y_r, omegaList,...
    gainList, S, ant_idx, numRefine)
% SUMMARY:
%   uses refineOne algorithm to refine frequencies and gains of
%   of all sinusoids
% INPUT:
%    y_r - residual measurement after all detected sinusoids have been
%          removed
%    omegaList - list of frequencies of known(detected) sinusoids
%    gainList  - list of gains of known(detected) sinusoids
%    S - measurement matrix
%        if [], measurements are direct i.e S = eye(N), where N is
%       length of sinusoid	
%    ant_idx - translating indexes to phases in definition of sinusoid
%    numRefine - number of times each sinusoid in the list is refined
%
% OUTPUT:
%       refined versions of inputs omegaList, gainList, y_r

K = length(omegaList); % number of sinusoids

% refinement repeated "numRefine" number of times per sinusoid

for i = 1:numRefine
    
    % chose an ordering for refinement
    
    % % *random* ordering
    % order = randperm(K);
    
    % *sequential* ordering
    order = 1:K;
    
    for j = 1:K
        l = order(j);
        
        % parameters of the l-th sinusoid
        omega = omegaList(l);
        gain = gainList(l);
        
        % refine our current estimates of (gain, omega) of the
        % l-th sinusoid
        [omega, gain, y_r] = refineOne(y_r,...
			       omega, gain, S, ant_idx, false);
        omegaList(l) = omega;
        gainList(l) = gain;
        % refineOne ensures that (gain, omega) pair are altered iff
        % the energy in the residual measurements y_r does not 
    	% increase
    end
    
end

end

% --------------------------------------------------------------------

function sampledManifold = preProcessMeasMat(S, overSamplingRate)
% SUMMARY:
%   compute overSamplingRate*M IFFTs once when measurement matrix
%   is not equal to eye(N); avoids repetitive IFFT computations 
%   later on
% INPUT:
%   S - M times N measurement matrix
%       M - number of measurements
%       N - length of sinusoid
%   examples: eye(N) - normal measurement matrix
%             diag(hamming(N)) - normal measurements with hamming
%                                weights
%             randn(M,N)/sqrt(N) - compressive measurement matrix
%             randn(M,N)/sqrt(N) * diag(hamming(N)) - compressive
%                       measurements with hamming weights matrix
%    overSamplingRate (optional) - how fine (in terms of multiples
%       of the FFT grid) do we want the coarse grid to be?
%       number of grid points = overSamplingRate * N
%       Default value is 3
%
% OUTPUT:
%       data structure with sinusoid responses 
%       precomputed using IFFTs when S is not equal to eye(N)


M = size(S,1);
N = size(S,2);

sampledManifold.length = N;
R = round(overSamplingRate*N);

sampledManifold.coarseOmega = 2*pi*(0:(R-1))/R;  % omegaCoarse

% definition of a sinusoid: 
%              exp(-1j*omega((N-1)/2):((N-1)/2))/sqrt(N)

ant_idx = 0:(N-1);
% ant_idx = ant_idx - (N-1)/2;

% IFFT definition of a sinusoid(omega) takes the following form:
% 	sinusoid    = @(omega) exp(1j*(0:(N-1)).'*omega);
% To reiterate, we assume that a sinusoid is given by
%	sinusoid    = @(omega) exp(1j*ant_idx.'*omega)/sqrt(N);
% So we store this information in sampledManifold container
sampledManifold.ant_idx = ant_idx(:);

% Check if the measurement matrix is the identity matrix
if M == N
    is_eye = norm(eye(N) - S, 'fro') == 0;
else
    is_eye = false;
end
sampledManifold.is_eye = is_eye;

% WHEN THE MEASUREMENT MATRIX NOT THE IDENTITY MATRIX
% we compute the Fourier transforms of sensing weights 
% and store them for use in the coarse stage (can use 
% IFFT to speed up these steps)

% IN THE FOLLOWING COMPUTATIONS WE COMPENSATE FOR THE 
% DIFFERENCES BETWEEN NOTION OF A SINUSOID THAT WE TAKE 
% AND THAT TAKEN BY "IFFT"

if ~is_eye
    
    % S times x(omegaCoarse)
    sampledManifold.map_IfftMat = R/sqrt(N)*ifft(S,R,2) * ...
        sparse(1:R, 1:R, exp(1j * sampledManifold.coarseOmega...
        * ant_idx(1)), R, R);
    
    % norm square of S times x(omegaCoarse)
    % energy as a function of frequency
    sampledManifold.map_IfftMat_norm_sq = ...
        sum(abs(sampledManifold.map_IfftMat).^2, 1);
    
end
end

 
