function [] = experiment()
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

 WITH_COMB  = false;
 ALGORITHM1 = true;
 VERBOSE    = false;
 % TIMING     = false;
 
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

  LARGE_FREQ = zeros(1,k);


%%%%%%%%% CALCULATING AND PLOTTING THE TRANSFORMS %%%%%%%%%  
%  for  idx=1:Dx:plotEnd 
    % for  idx=1:Dx:plotEnd/3 

  for  idx=1:Dx:Dx 
      clf;
      
      R1=real_1(idx : idx+Dx-1);
      I1 = imag_1(idx : idx+Dx-1);
      x = R1 + 1i*I1;

      sfft_real = fftshift(run_experiment(R1, n,
                    lobefrac_loc, tolerance_loc, b_loc,
                    B_loc, B_thresh, loc_loops, threshold_loops,
                    lobefrac_est, tolerance_est, b_est,
                    B_est, est_loops, W_Comb, Comb_loops,
                    repetitions, FFTW_OPT, LARGE_FREQ, k));

      sfft_imag = fftshift(run_experiment(I1, n,
                    lobefrac_loc, tolerance_loc, b_loc,
                    B_loc, B_thresh, loc_loops, threshold_loops,
                    lobefrac_est, tolerance_est, b_est,
                    B_est, est_loops, W_Comb, Comb_loops,
                    repetitions, FFTW_OPT, LARGE_FREQ, k));
      sfft = fftshift(run_experiment(x, n,
                    lobefrac_loc, tolerance_loc, b_loc,
                    B_loc, B_thresh, loc_loops, threshold_loops,
                    lobefrac_est, tolerance_est, b_est,
                    B_est, est_loops, W_Comb, Comb_loops,
                    repetitions, FFTW_OPT, LARGE_FREQ, k));

      
      dB_cmplx = 20*log10(3.4*abs(sfft));  
      dB_real = 20*log10(abs(sfft_real));
      dB_imag = 20*log10(abs(sfft_imag));
      dBFS_real_s = dB_real - 20*log10(256)-20*log10(2^15)+20*log10(2^(0.5))-2.0;
      dBFS_imag_s = dB_imag - 20*log10(256)-20*log10(2^15)+20*log10(2^(0.5))-2.0;
      dBFS_cmplx_s = dB_cmplx - 20*log10(256)-20*log10(2^15)+20*log10(2^(0.5))-2.0;
      
      subplot(2,2,1)
      hold all;
      plot( distance,R1,'linewidth',1);
      plot( distance,I1,'linewidth',1);
      title('IQ time domain signal','FontSize',20);
      % axis([-25 25 -500 500]);
      xlabel('Samples','FontSize',18);
      ylabel('ADC value','FontSize',18);
      grid on ;
      set(gca,'FontSize',18,'FontWeight','bold');

      subplot(2,2,2);  
      grid on ;
      hold all;
      stem(distance,dBFS_real_s,'linewidth',2);
      %    axis([-25 25 -140 0]);
      title('SFFT Amplitude(per chirp-real only)','FontSize',20);
      xlabel('Distance[m]','FontSize',18);
      ylabel('FFT output [dBFS]','FontSize',18);
      set(gca,'FontSize',18,'FontWeight','bold');


      subplot(2,2,3);
      grid on ;
      hold all;
      stem(distance,dBFS_imag_s,'linewidth',2);
      % axis([-25 25 -140 0]);
      title('SFFT Amplitude(per chirp-imag only)','FontSize',20);
      xlabel('Distance[m]','FontSize',18);
      ylabel('FFT output [dBFS]','FontSize',18);
      set(gca,'FontSize',18,'FontWeight','bold');

      subplot(2,2,4);
      grid on ;
      hold all;
      stem(distance,dBFS_cmplx_s,'linewidth',1);
      % axis([-25 25 -140 0]);
      title('SFFT Amplitude(per chirp)','FontSize',20);
      xlabel('Distance[m]','FontSize',18);
      ylabel('FFT output [dBFS]','FontSize',18);
      set(gca,'FontSize',18,'FontWeight','bold');

      R11 = fftshift(fft_recur(hanning(length(R1)).*R1));
      I11 = fftshift(fft_recur(hanning(length(I1)).*I1));
      ABS11 = fftshift(fft_recur(hanning(length(R1+1i*I1)).*(R1+1i*I1)));
      dB_cmplx = 20*log10(abs(ABS11));  
      dB_real = 20*log10(abs(R11));
      dB_imag = 20*log10(abs(I11));
      dBFS_real = dB_real  - 20*log10(256)-20*log10(2^15)+20*log10(2^(0.5))-2.0;
      dBFS_imag = dB_imag - 20*log10(256)-20*log10(2^15)+20*log10(2^(0.5))-2.0;
      dBFS_cmplx= dB_cmplx - 20*log10(256)-20*log10(2^15)+20*log10(2^(0.5))-2.0;

      % figure;
      subplot(2,2,1)
      hold all;
      plot( distance,R1,'b','linewidth',1);
      plot( distance,I1,'r','linewidth',1);
      title('IQ time domain signal','FontSize',20);
      axis([-25 25 -500 500]);
      xlabel('Samples','FontSize',18);
      ylabel('ADC value','FontSize',18);
      grid on ;
      set(gca,'FontSize',18,'FontWeight','bold');
    
      subplot(2,2,2);  
      grid on ;
      hold all;
      plot(distance,dBFS_real,'k','linewidth',2);
      % axis([-25 25 -140 0]);
      title('1D FFT Amplitude profile(per chirp-real only)','FontSize',20);
      xlabel('Distance[m]','FontSize',18);
      ylabel('FFT output [dBFS]','FontSize',18);
      set(gca,'FontSize',18,'FontWeight','bold');
      legend('SFFT','Normal FFT','FontSize',7)

      
      subplot(2,2,3);
      grid on ;
      hold all;
      plot(distance,dBFS_imag,'g','linewidth',2);
      %  axis([-25 25 -140 0]);
      title('1D FFT Amplitude profile(per chirp-imag only)','FontSize',20);
      xlabel('Distance[m]','FontSize',18);
      ylabel('FFT output [dBFS]','FontSize',18);
      set(gca,'FontSize',18,'FontWeight','bold');
      legend('SFFT','Normal FFT','FontSize',7)

      subplot(2,2,4);
      grid on ;
      hold all;
      plot(distance,dBFS_cmplx,'r','linewidth',1);
      %  axis([-25 25 -140 0]);
      title('1D FFT Amplitude profile(per chirp)','FontSize',20);
      xlabel('Distance[m]','FontSize',18);
      ylabel('FFT output [dBFS]','FontSize',18);
      set(gca,'FontSize',18,'FontWeight','bold');
      legend('SFFT','Normal FFT','FontSize',7)
      drawnow;
  end
end

%%%%%%%%% CALCULATING FILTER VALUES AND RUN OUTER LOOP  %%%%%%%%%  
function [x_f] = run_experiment(x, n, lobefrac_loc, tolerance_loc, b_loc, B_loc, B_thresh, loops_loc, loops_thresh,
		    lobefrac_est, tolerance_est, b_est, B_est, loops_est, W_Comb, Comb_loops,
		    repetitions, bool FFTW_OPT, LARGE_FREQ, k)
  
  [filtert,w_loc] = make_dolphchebyshev_t(lobefrac_loc, tolerance_loc);
  [filter_timedo,filter_sizet,filter_freqdo] = make_multiple_t(filtert, w_loc, n, b_loc);
  

  [filtert_est,w_est] = make_dolphchebyshev_t(lobefrac_est, tolerance_est, w_est);
  [filter_est_timedo,filter_est_sizet,filter_est_freqdo] = make_multiple_t(filtert_est, w_est, n, b_est);
  
  x_f = outer_loop(x, n, filter_timedo,filter_sizet,filter_freqdo, filter_est_timedo,
            filter_est_sizet,filter_est_freqdo, B_est, B_thresh, B_loc, W_Comb,
            Comb_loops, loops_thresh, loops_loc, loops_loc + loops_est);
        
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
    h(mod((ii+n +offset),n)+1) = s;
    maximum = max(maximum,abs(s));
    s = s + (g(mod((ii + b),n)+1) - g(ii+1));
  end

  h = h/maximum;

  offsetc = 1;
  const_gain=exp((-2*pi * 1i * (w/2)) / n);
  
  for ii = 0:n-1
    h(ii+1) = h(ii+1)*offsetc;
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

function answer = timesmod(x, a, n) 
   answer = int64(mod((int64(x) * a),n));
end

function v = mod_inverse(a, n) 
 ii = n;
 v = 0;
 d = 1;
 
 while a>0 
  t = ii/a;
  x = a;
  
  a = mod(ii , x);
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
