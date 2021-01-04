function [] = experiment()
 
 n = 256;
 k = 17;
 repetitions = 1;
 Bcst_loc = 32;
 Bcst_est = 32;
 Comb_cst = 16;
 loc_loops = 0;
 est_loops = 6;
%  threshold_loops = round((loc_loops+est_loops)/2);
 threshold_loops = 3;
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
  lobefrac_est = 0.2 / (BB_est);

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
  fprintf("  threshold_loops : %d \n ",threshold_loops);
  fprintf("  B_loc : %d \n ",B_loc);
  fprintf("  B_thresh : %d \n ",B_thresh);
  fprintf("  B_est : %d \n ",B_est);

  x  = rand(1,n);

    % disp("SFFT REAL BEING CALCULATED");
    sfft = fftshift(run_experiment(x',256,lobefrac_loc, tolerance_loc, b_loc,B_loc, B_thresh, loc_loops, threshold_loops,lobefrac_est, tolerance_est, b_est,B_est, est_loops, W_Comb, Comb_loops,repetitions, FFTW_OPT, LARGE_FREQ, k));

    [sorted,I] = sort(abs(sfft),'descend');
    % I = (abs(sfft)~=0);
    temp = zeros(1,length(sfft));
    % temp(I) = ABS11(I);
    % temp(I(1:k)) = ABS11(I(1:k));
    temp(I(1:k)) = sfft(I(1:k));
    sfft = temp;
    
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
  
  w = (1 / pi) * (1/lobefrac) * acosh(1./tolerance);
  % disp(w);
  w=255;
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
  const_gain=exp((-2*pi * 1i * double(w/2)) / double(n));
  % const_gain=1;
  
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
