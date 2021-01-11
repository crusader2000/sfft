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
fid = fopen('../adc_data_Raw.bin','r');
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
 

%  figure(1);
%  hold on;

 WITH_COMB  = false;
 ALGORITHM1 = true;
 VERBOSE    = false;
 % TIMING     = false;
 
 n = 256;
 k = 1;
 repetitions = 1;
 Bcst_loc = 32;
 Bcst_est = 32;
 Comb_cst = 16;
 loc_loops = 0;
 est_loops = 3;
%  threshold_loops = round((loc_loops+est_loops)/2);
 threshold_loops = 1;
 Comb_loops = 1;
 simulate = 0;
 snr = 1000000000;
 std_noise = 0;
 FFTW_OPT = false;
 tolerance_loc = 1.e-6;
 tolerance_est = 1.e-8; 
 n = floor_to_pow2(n);

  BB_loc =  uint8(Bcst_loc*sqrt(n*k/(log2(n))));
  BB_est =  uint8(Bcst_est*sqrt(n*k/(log2(n))));

  lobefrac_loc = 0.5 / (BB_loc);
  lobefrac_est = 0.5 / (BB_est);

  b_loc = int32(1.2*1.1*( n/BB_loc));
  b_est = int32(1.4*1.1*( n/BB_est));

  B_loc = floor_to_pow2(BB_loc);
  B_thresh = 2*k;
  B_est = floor_to_pow2(BB_est);

  W_Comb = floor_to_pow2(Comb_cst*n/B_loc);

  LARGE_FREQ = zeros(1,k);

  fprintf("  k : %d \n ",k);
  fprintf("  loc_loops : %d \n ",loc_loops);
  fprintf("  est_loops : %d \n ",est_loops);
  fprintf("  est_loops : %d \n ",est_loops);
  fprintf("  threshold_loops : %d \n ",threshold_loops);
  fprintf("  B_loc : %d \n ",B_loc);
  fprintf("  B_thresh : %d \n ",B_thresh);
  fprintf("  B_est : %d \n ",B_est);

  % disp(B_loc);
  % disp(B_thresh);
  % disp(B_est);

  runtime_fft = 0;
  runtime_sfft = 0;

  chirp_count = 0;
  abs_acc = 0;
  abs_acc_arr = [];
  abs_acc_phase_arr = [];
n = 256;
  tic;

%%%%%%%%% CALCULATING AND PLOTTING THE TRANSFORMS %%%%%%%%%  
%  for  idx=1:Dx:plotEnd 
% for est_loops = 3:3

  chirp_count = 0;
  abs_acc = 0;
  abs_phase_acc = 0;
for  idx=1:Dx:plotEnd/3
  % for  idx=1:Dx:Dx  
    % for  idx=Dx+1:Dx:2*Dx  

    chirp_count = chirp_count + 1;

    % w = waitforbuttonpress;

    clf;
    R1=real_1(idx : idx+Dx-1);
    I1 = imag_1(idx : idx+Dx-1);
    x = R1 + 1i*I1;
    
    time_temp = toc;
    ABS11 = fftshift(fft_recur(hanning(length(R1+1i*I1)).*(R1+1i*I1)));
    runtime_fft = toc - time_temp;
   
    % dB_cmplx = 20*log10(abs(ABS11));  
    
    dB_cmplx = abs(ABS11);  
    
    % dB_cmplx = angle(ABS11);  

    % dBFS_cmplx = dB_cmplx - 20*log10(256)-20*log10(2^15)+20*log10(2^(0.5))-2.0;
    dBFS_cmplx = dB_cmplx;


    % disp("SFFT REAL BEING CALCULATED");
    time_temp = toc;
    sfft = fftshift(run_experiment(x',256,lobefrac_loc, tolerance_loc, b_loc,B_loc, B_thresh, loc_loops, threshold_loops,lobefrac_est, tolerance_est, b_est,B_est, est_loops, W_Comb, Comb_loops,repetitions, FFTW_OPT, LARGE_FREQ, k));


    disp("ABS VALUES");
    [M,I] = max(abs(sfft));
    temp = zeros(1,length(sfft));
    temp(I) = sfft(I);
    A = I;
    % A_phase = angle(temp(I));
    [M,I] = max(abs(ABS11));
    B = I;
    factor_dista =((10^7)/256)*((1.5*10^8)/(29.982*10^12));

    fprintf("A = %d   B = %d\n",A,B);
    fprintf("A = %d   B = %d (in metres)\n",(A-128)*factor_dista,(B-128)*factor_dista);
    fprintf("f(A) = %d   f(B) = %d\n",abs(sfft(A)),abs(ABS11(B)));

    % B_phase = angle(ABS11(I));
    
    % abs_acc = abs_acc + get_freq_mse(A,B);
    % abs_phase_acc = abs_phase_acc + get_phase_mse(A,A_phase,B,B_phase);
    sfft = temp;
    
    runtime_sfft = toc - time_temp;
    
    % dB_cmplx = 20*log10(abs(sfft));  
    
    dB_cmplx = abs(sfft);  
    
    % % PHASE OF THE SIGNAL
    % dB_cmplx = angle(sfft);  
    
    % dBFS_cmplx_s = dB_cmplx - 20*log10(256)-20*log10(2^15)+20*log10(2^(0.5))-2.0;

    dBFS_cmplx_s = dB_cmplx;

    grid on ;
    hold all;
    stem(distance,dBFS_cmplx_s,'linewidth',1);
    % axis([-25 25 -140 0]);
    title('SFFT Amplitude(per chirp)','FontSize',20);
    xlabel('Distance[m]','FontSize',18);
    ylabel('FFT output','FontSize',18);
    set(gca,'FontSize',18,'FontWeight','bold');

    % figure;

    grid on ;
    hold all;
    plot(distance,dBFS_cmplx,'r','linewidth',1);
    %  axis([-25 25 -140 0]);
    title('1D FFT Amplitude profile(per chirp)','FontSize',20);
    xlabel('Distance[m]','FontSize',18);
    ylabel('FFT output','FontSize',18);
    set(gca,'FontSize',18,'FontWeight','bold');
    legend('SFFT','Normal FFT','FontSize',7)
    drawnow;
    % figure;

    % grid on ;
    % hold all;
    % stem(distance,angle(sfft),'linewidth',1);
    % % axis([-25 25 -140 0]);
    % title('SFFT Amplitude(per chirp)','FontSize',20);
    % xlabel('Distance[m]','FontSize',18);
    % ylabel('FFT output','FontSize',18);
    % set(gca,'FontSize',18,'FontWeight','bold');

    % % figure;

    % grid on ;
    % hold all;
    % plot(distance,angle(ABS11),'r','linewidth',1);
    % %  axis([-25 25 -140 0]);
    % title('1D FFT Amplitude profile(per chirp)','FontSize',20);
    % xlabel('Distance[m]','FontSize',18);
    % ylabel('FFT output','FontSize',18);
    % set(gca,'FontSize',18,'FontWeight','bold');
    % legend('SFFT','Normal FFT','FontSize',7)
    % drawnow;

    % figure;

% pause(2);
  end

  % abs_acc_arr(est_loops-1) =  sqrt(abs_acc/chirp_count);
  % abs_acc_phase_arr(est_loops-1) = sqrt(abs_phase_acc/chirp_count);

  % fprintf("RMSE of ABS location/frequency %d \n", sqrt(abs_acc/chirp_count));
  % fprintf("RMSE of ABS phases (in radians) %d \n", sqrt(abs_acc/chirp_count));
  % fprintf("RMSE of ABS phases (in degrees) %d \n", sqrt(abs_acc/chirp_count)*(180/3.14));  

% end
  % figure;
  % plot(1:length(abs_acc_arr),abs_acc_arr);
  % title("Frequency RMSE (Imag)");
  % xlabel("No. of outer loops");
  % ylabel("RMSE");
  % figure;
  % plot(1:length(abs_acc_phase_arr),abs_acc_phase_arr);
  % title("Phase RMSE (Imag)");
  % xlabel("No. of outer loops");
  % ylabel("RMSE (Radians)");
  
  % fprintf("Average accuracy of ABS values %d \n", (abs_acc/chirp_count)*100);

  % fprintf(" Runtime of FFT is %d \n ",runtime_fft);
  % fprintf(" Runtime of SFFT is %d \n \n",runtime_sfft);
end

function acc = get_freq_mse(A,B)
  num1 = length(A);
  num2 = length(B);
  acc = 0;
  for ii = 1:num1
    for jj = 1:num2
      if A(ii) == B(jj)
        % acc = acc + 1;
        A(ii) = 0;
        B(jj) = 0;
        break; 
      end
    end
  end
 
  A = sort(A);
  B = sort(B);
  % disp("A");
  % disp(A);
  % disp("B");
  % disp(B);
  acc = mean((A'-B).^2);
  % disp(acc);
  % acc = acc/num1;
  % if acc == 0
  % disp("Accuracy Zero");
  % disp(A);
  % disp(B);
  % end
end

function acc = get_phase_mse(A,A_phase,B,B_phase)
  num1 = length(A);
  num2 = length(B);
  acc = 0;
  
  for ii = 1:num1
    for jj = 1:num2
      if A(ii) == B(jj)
        % acc = acc + 1;
        A(ii) = 0;
        B(jj) = 0;
        acc = acc + (A_phase(ii)-B_phase(jj))^2;
        break; 
      end
    end
  end
 
  [A,I] = sort(A);
  B = sort(B);

  for ii = 1:num1
    if A(ii) ~= 0
      acc = acc + (A_phase(I(ii))-B_phase(I(ii)))^2;
    end
  end
  % disp("A");
  % disp(A);
  % disp("B");
  % disp(B);
  % acc = mean((A'-B).^2);
  % disp(acc);
  acc = acc/num1;
  % if acc == 0
  % disp("Accuracy Zero");
  % disp(A);
  % disp(B);
  % end
end

%%%%%%%%% CALCULATING FILTER VALUES AND RUN OUTER LOOP  %%%%%%%%%  
function x_f = run_experiment(x, n, lobefrac_loc, tolerance_loc, b_loc, B_loc, B_thresh, loops_loc, loops_thresh,lobefrac_est, tolerance_est, b_est, B_est, loops_est, W_Comb, Comb_loops,repetitions,FFTW_OPT, LARGE_FREQ, k)
  
  [filtert,w_loc] = make_dolphchebyshev_t(lobefrac_loc, tolerance_loc);
  [filter_timedo,filter_sizet,filter_freqdo] = make_multiple_t(filtert, w_loc, n, b_loc);

  [filtert_est,w_est] = make_dolphchebyshev_t(lobefrac_est, tolerance_est);
  [filter_est_timedo,filter_est_sizet,filter_est_freqdo] = make_multiple_t(filtert_est, w_est, n, b_est);

  fprintf("Filter length  w = %d \n",filter_sizet);
  fprintf("Filter Est length  w_est = %d \n",filter_est_sizet);

  % plot(1:length(filter_timedo),abs(filter_timedo));
  % figure;

  % disp(length(filter_freqdo));
  % disp("CALCULATED FILTERS");  

  x_f = outer_loop(x, n, filter_timedo,filter_sizet,filter_freqdo, filter_est_timedo,filter_est_sizet,filter_est_freqdo, B_est, B_thresh, B_loc, W_Comb,Comb_loops, loops_thresh, loops_loc, loops_loc + loops_est);
  
  % plot(1:length(filter_freqdo),20*log10(abs(filter_freqdo)));
  % figure;
end

% WINDOW FUNCTIONS
function [x,w] = make_dolphchebyshev_t(lobefrac,tolerance)
  % disp("CALCULATING DOLPHCHEBYSHEV FILTER");  
  
  w = double(((1 / pi) * (1/lobefrac) * acosh(1./tolerance)));
  % w=255;
  disp("w");
  disp(w);
  disp(class(w));
  if ~mod(w,2) 
   w = w+1;
 end

  x = zeros(1,w);

  t0 = cosh(double(acosh(1/tolerance) / (w-1)));

 
 for ii=0:w-1
    x(ii+1) = Cheb(w-1, t0 * cos(double((pi * ii) / w))) * tolerance;
  end

  x = fft(x);
  %  x = fftshift(x);

 x= real(x);
end

function out = Cheb(m, x)
  if abs(x) <= 1
    out = cos(double(m * acos(x)));
  else
    out = real(cosh(m * acosh(x)));
  end
end

function [x,w,h] = make_multiple_t(x,w,n,b)
  % x - Time Domain
  % w - Size
  % h - Frequency Domain
  g = zeros(1,n);
  h = zeros(1,n);
 
 
  g(1:round(w/2)) = x(round(w/2):w);
  g((n - round(w/2)):uint8(n)) = x(1:round(w/2));

  g = fft(g);
  
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
  % const_gain=exp((-2*pi * 1i * double(w/2)) / double(n));
  const_gain=1;
  
  for ii = 0:n-1
    h(ii+1) = h(ii+1)*offsetc;
    offsetc = offsetc*const_gain;
  end

  % g = fft(h);
  % x = g(1:w)/n;
end


% UTILITY FUNCTIONS
function answer = floor_to_pow2(x)
    ii = 1;
    while ii<= x
      answer = ii;
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
   answer = int32(mod((int32(x) * a),n));
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
