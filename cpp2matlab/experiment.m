function x_f = experiment()
% x_f - True FFT
clc;
clear all;
close all;
%% global variables
% Based on sensor configuration.
   numADCBits = 16; % number of ADC bits per sample.
   numADCSamples = 256; % number of ADC samples per chirp.
   numRx = 4; % number of receivers in AWR1243.
   chirpSize = numADCSamples*numRx;
   chirploops= 128; % No. of of chirp loops.  
   numLanes = 4; % do not change. number of lanes is always 4
   isReal = 0; % set to 1 if real only data, 0 if complex data.
   numFrames = 200; 
   numChirps = 1;% To consider all the chriploops keep numchrirps = chirploops.
   sampleRate = 10; % [Msps]
   timeStep = 1/sampleRate;    % [us]
   chirpPeriod = numADCSamples * timeStep ; % [us]
   plotEnd = numADCSamples *  numChirps *numFrames; %for considering all frames.
   Dx = numADCSamples * numChirps ;
   timeEnd = (plotEnd-1) * timeStep;

%% read file
% read .bin file
fid = fopen('../adc_data_Raw_Raw_0.bin','r');
adcData = fread(fid, 'int16');
fclose(fid);
% adcData = adcData1(1:length(adcData1)/numFrames);
fileSize = size(adcData, 1);
% one chirp only (next line iin case of multiple chirps)
% adcData = adcData(1:chirpSize,1);
% organize data by LVDS lane
% for complex data
  remaind = mod(fileSize,8);
% Make data(Interleaved Data from AWR1243) over 8 columns. 
if remaind ~= 0
   adcData =[ adcData;zeros(8-remaind,1)] ;
end
   fileSize = length(adcData);
%% stroing data in LVDS if Real and in cmplx if complex(IQ from mmwave studio)   
if isReal % For real data 4 columns for 4 receivers
    adcData = adcData';
    LVDS = reshape(adcData ,4,[])';

else
% cmplx has 4 real & 4 imaginary columns for 4 Rceivers for interleaved data format.
    adcData = adcData';
    cmplx = reshape(adcData ,8,[])';
end

%% return receiver data
if isReal 
    retValue = LVDS;
else
    retValue = cmplx;
end

%% plotting the data
adcData = retValue ;

sample = (0:1:plotEnd-1);
time = (0:timeStep:timeEnd);
f = (0:1:plotEnd-1);
f_bin = (0:1:length(adcData)-1);

% % Distance calculation using d=(c*f/(2*slope))
fdel_bin = (-128:1:127)*((10^7)/256);
distance = ((1.5*10^8)*fdel_bin)/(29.982*10^12);


% % plotting for all frames (Channel 1)
 
 real_1 = adcData(:,1);        
 imag_1 = adcData(:,5);
 
 figure(1);
 hold on;

n = 256;
k = 20;
repetitions = 1;
Bcst_loc=2;
Bcst_est=0.2;
Comb_cst=16;
loc_loops =5;
est_loops =12;
threshold_loops =3;
Comb_loops = 1;
simulate = 0;
snr=1000000000;
std_noise = 0;
FFTW_OPT = false;
tolerance_loc = 1.e-6;
tolerance_est = 1.e-8; 
n = floor_to_pow2(n);

BB_loc =  uint8(Bcst_loc*sqrt(n*k/(log2(n))));
BB_est =  uint8(Bcst_est*sqrt(n*k/(log2(n))));

lobefrac_loc = 0.5 / (BB_loc);
lobefrac_est = 0.5 / (BB_est);

b_loc = int64(1.2*1.1*( n/BB_loc));
b_est = int64(1.4*1.1*( n/BB_est));

B_loc = floor_to_pow2(BB_loc);
B_thresh = 2*k;
B_est = floor_to_pow2(BB_est);

W_Comb = floor_to_pow2(Comb_cst*n/B_loc);

x = zeros(1,n);

x_f = zeros(1,n);

LARGE_FREQ = zeros(1,k);

x = fft_recur(x_f);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
INPUT DATA ISSUES
x_f = run_experiment(x, n,
               lobefrac_loc, tolerance_loc, b_loc,
               B_loc, B_thresh, loc_loops, threshold_loops,
               lobefrac_est, tolerance_est, b_est,
               B_est, est_loops, W_Comb, Comb_loops,
               repetitions, FFTW_OPT, LARGE_FREQ, k);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

function [x_f] = run_experiment(x, n, lobefrac_loc, tolerance_loc, b_loc, B_loc, B_thresh, loops_loc, loops_thresh,
		    lobefrac_est, tolerance_est, b_est, B_est, loops_est, W_Comb, Comb_loops,
		    repetitions, bool FFTW_OPT, LARGE_FREQ, k)
  
  
  [filtert,w_loc] = make_dolphchebyshev_t(lobefrac_loc, tolerance_loc);
  [filter_timedo,filter_sizet,filter_freqdo] = make_multiple_t(filtert, w_loc, n, b_loc);
  

  [filtert_est,w_est] = make_dolphchebyshev_t(lobefrac_est, tolerance_est, w_est);
  [filter_est_timedo,filter_est_sizet,filter_est_freqdo] = make_multiple_t(filtert_est, w_est, n, b_est);
  

  % filter_noise = 0;
  % filter_noise_est = 0;
  
  % for(int i = 0; i < 10; i++) {
  %   filter_noise = std::max(filter_noise,
	% 		    std::max(cabs(filter.freq[n/2+i]),
	% 			     cabs(filter.freq[n/2-i])));
  %   filter_noise_est = std::max(filter_noise_est,
  %                               std::max(cabs(filter_est.freq[n/2+i]),
  %                                        cabs(filter_est.freq[n/2-i])));
  % }
  % printf(" Noise in filter: Location Filter : %lg; Estimation Filter %lg\n", filter_noise, filter_noise_est);
  % printf("******************************************************************************\n\n");


  % printf("sFFT Results\n");
  % printf("******************************************************************************\n");

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
INPUT DATA ISSUES
x_f = outer_loop(origx,n, filter_timedo,filter_sizet,filter_freqdo, 
                        filter_est_timedo,filter_est_sizet,filter_est_freqdo, B2,
                        num,B,W_Comb,Comb_loops,loop_threshold,location_loops,loops)



  % int num_candidates= (int)oloop_output.size();
  % std::pair<real_t, int> *candidates = (std::pair<real_t, int> *)malloc(num_candidates*sizeof(*candidates));
  % complex_t *x_f_Large = (complex_t *)calloc(n,sizeof(*x_f_Large));
  % complex_t *ans_Large = (complex_t *)calloc(n,sizeof(*ans_Large));

  % int counter=0;

  % for(__typeof(oloop_output.begin()) it = oloop_output.begin(); it != oloop_output.end(); it++){
  %    int key = it->first;
  %    complex_t value =it->second;
  %    candidates[counter] = std::make_pair(cabs(value), key);
  %    counter++;
  %  }

  % //Enter ALL large frequences as zero
  % for(int i = 0; i < k; i++){
  %     x_f_Large[LARGE_FREQ[i]]=x_f[LARGE_FREQ[i]];
  % }

  % std::nth_element(candidates, candidates + num_candidates - k, candidates + num_candidates);
  % for(int i = 0; i < k; i++){
  %   int key = candidates[num_candidates - k + i].second;
  %   ans_Large[key] = oloop_output[key];
  % }

  % int large_found = 0;
  % int FOUND=0;
  % for(int i = 0; i < k; i++){
  %   FOUND += (unsigned int)oloop_output.count(LARGE_FREQ[i]);
  %   large_found += (ans_Large[(LARGE_FREQ[i])] != 0);
  % }

end

% WINDOW FUNCTIONS
function [x,w] = make_dolphchebyshev_t(lobefrac,tolerance)
  w = int((1 / pi) * (1/lobefrac) * acosh(1./tolerance));
 if ~(mod(w,2))
   w = w-1;
 end

 x = zeros(1,w);
 t0 = cosh(acosh(1/tolerance) / (w-1));
 
 for ii=0:w-1
  x(ii+1) = Cheb(w-1, t0 * cos((pi * ii) / w)) * tolerance;
  end
 x = fft_recur(x);
 fftshift(x);

 x= real(x);
end

function out = Cheb(m, x)
  if abs(x) <= 1
    out = cos(m * acos(x));
  else
    out = real(ccosh(m * cacosh(x)));
  end
end

function [x,w,h] = make_multiple_t(x,w,n,b)
  % x - Time Domain
  % w - Size
  % h - Frequency Domain

  g = zeros(1,n);
  h = zeros(1,n);
 
 
  % memcpy(g, x+w/2, (w - (w/2))*sizeof(*g));
  % memcpy(g + n - w/2, x, (w/2)*sizeof(*g)); 
  g(1:(w/2)) = x(w/2+1:w);
  g((n - (w/2) +1):n) = x(1:w/2);


  g = fft_recur(g);
  
  s = sum(g(1:b));


  maximum = 0;
  offset = b/2;

  for ii = 0:n-1
    h(mod((i+n +offset),n)+1) = s;
    maximum = max(maximum,abs(s));
    s = s + (g(mod((i + b),n)+1) - g(i+1));
  end

  h = h/maximum;

  offsetc = 1;
  const_gain=exp((-2*pi * I * (w/2)) / n);
  
  for ii = 0:n-1
    h(i) = h(i)*offsetc;
    offsetc = offsetc*const_gain;
  end


  g = fft_recur(h);
  x = g(1:w)/n;


end


% UTILITY FUNCTIONS
function answer = floor_to_pow2(x)
    for ii = 1:x
        if ii<= x
            answer = ii;
        end
        ii = 2*ii;
    end
    answer = answer/2;
end

function answer = gcd(a,b)
  if (mod(a,b) == 0) 
    answer =  b;
  else 
    answer = gcd(b, mod(a,b)); 
  end
end

function answer = timesmod(const int &x, const int &a, const int &n) 
   answer = int64(mod((int64(x) .* a),n));
end

function v = mod_inverse(a, n) 
 ii = n;
 v = 0;
 d = 1;
 
 while a>0 
  t = ii/a;
  x = a;
  
  a = ii % x;
  ii = x;
  x = d;
  d = v - t*x;
  v = x;
 end

 v = mod(v,n);
 if v<0 
  v = mod((v+n),n);
 end
end

function output = find_largest_indices(num,samples)
 [sorted,I] = sort(samples);
 output = I(1:num);
end 
