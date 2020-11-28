function oloop_output = outer_loop(origx,n, filter_timedo,filter_sizet,filter_freqdo, 
                                  filter_est_timedo,filter_est_sizet,filter_est_freqdo, B2,
                                  num,B,W_Comb,Comb_loops,loop_threshold,location_loops,
                                  loops)

    permute_vals = zeros(1,loops);
    permuteb_vals = zeros(1,loops);

  
  x_samp = {}
  for ii = 0:loops
    if ii<location_loops
    x_samp{ii+1} = zeros(1,B);
    else
    x_samp{ii+1} = zeros(1,B2);
    end
  end
 
  hits_found = 0;

  if WITH_COMB
    pages_hit = num * (n/B) * num * (1 / W_Comb) * location_loops;
  else
    pages_hit = num * (n/B) * location_loops;
  end

 score = zeros(1,n);
 hits = zeros(1,n);
  
  %  BEGIN INNER LOOPS

  for ii=0:loops-1
    a = 0;
    b = 0;

    while gcd(a, n) ~= 1
      a = int(mod(randi(n),n));
    end

    ai = mod_inverse(a, n);

    permute_vals(ii+1) = ai;
    permuteb_vals(ii+1) = b;

    perform_location = int64(ii < location_loops);

    if perform_location
      cur_B = B;
      cur_filter_timedo = filter_timedo;
      cur_filter_sizet = filter_sizet;
      cur_filter_freqdo = filter_freqdo;
    else
      cur_B = B2;
      cur_filter_timedo = filter_est_timedo;
      cur_filter_sizet = filter_est_sizet;
      cur_filter_freqdo = filter_est_freqdo;  
    end

    J = zeros(1,num);

    [x_samp{ii+1},J] = inner_loop_locate(origx, n, cur_filter_timedo,
                      cur_filter_sizet,cur_filter_freqdo, num, cur_B,
                      a, ai, b, J);

    if perform_location 
      [score, hits, hits_found] = inner_loop_filter_regular(J, n, num, cur_B,
                                  a, ai, b, loop_threshold,
                                  score, hits, hits_found);
    end

  
  end
  
  %END INNER LOOPS
  oloop_output = estimate_values(hits,hits_found,
      x_samp, loops,n, permute_vals, B, B2,
      filter_timedo,filter_sizet,filter_freqdo, 
      filter_est_timedo,filter_est_sizet,filter_est_freqdo,location_loops)
end


function [x_samp,J] = inner_loop_locate(origx,n, filter_timedo,filter_sizet,filter_freqdo,
                     num,B,a,ai,b)

  if mod(n,B) ~= 0
    disp("Warning: n is not divisible by B, which algorithm expects.\n");
  end
  
  x_sampt = zeros(1,n);

  index=b;
  for ii = 0:filter_sizet
    x_sampt[mod(ii,B)+1] += origx[index+1] * filter_timedo[ii+1];
    index = mod((index+ai),n);
  end

  x_samp = fft_recur(x_sampt);

  samples = abs(x_samp[:B]);

  J= find_largest_indices(num,samples);
end


function [score, hits, hits_found] = inner_loop_filter_regular(J,n,num,B,a,ai,b,loop_threshold,
                              score, hits, hits_found)
  for ii=0:num-1
    low = mod((int64(ceil((J(ii+1) - 0.5) * n / B)) + n),n);
    high = mod((int64(ceil((J(ii+1) + 0.5) * n / B)) + n),n);
    loc = timesmod(low, a, n);
    
    jj= low;
    while jj ~= high
      score(loc+1) = score(loc+1)+1;

        if score(loc+1)==loop_threshold 
          hits(hits_found+1)=loc;
          hits_found = hits_found + 1;
          loc = mod((loc + a),n);
        end

      jj = mod((jj+1),n);
    end
  end
end

function oloop_output = estimate_values(hits,hits_found,
                x_samp, loops,n, permute_vals, B, B2,
                filter_timedo,filter_sizet,filter_freqdo, 
                filter_est_timedo,filter_est_sizet,filter_est_freqdo,location_loops)
  
  oloop_output = zeros(1,n);
  values = zeros(2,loops);

  for ii = 0:hits_found-1
    position = 0;

    for jj = 0:loops-1
      if jj <location_loops
        cur_B = B;
        cur_filter_timedo = filter_timedo;
        cur_filter_sizet = filter_sizet;
        cur_filter_freqdo = filter_freqdo;
      else
        cur_B = B2;
        cur_filter_timedo = filter_est_timedo;
        cur_filter_sizet = filter_est_sizet;
        cur_filter_freqdo = filter_est_freqdo;  
      end
      permuted_index= timesmod(permute_vals(jj+1), hits(ii+1),  n);
      hashed_to = (permuted_index / (n / cur_B));
      dista = mod(permuted_index , (n / cur_B));
     
      if dista > ((n/cur_B)/2) 
        hashed_to = mod((hashed_to + 1),cur_B);
        dista = (dista -(n/cur_B));
      end
      dista = mod((n - dista),n);

      filter_value = cur_filter_freqdo(dista);
      
      values(1,position+1) = real(x_samp{jj+1}(hashed_to) / filter_value);
      values(2,position+1) = imag(x_samp{jj+1}(hashed_to) / filter_value);
      position = position +1;

    end

    location = int64((loops + 1) / 2);

    values(1,:) = median(values(1,:));
    values(2,:) = median(values(2,:));
    
    realv = values(1,location);
    imagv = values(2,location);

    oloop_output(hits(ii+1)) = realv + 1i*imagv;
    
  end
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
