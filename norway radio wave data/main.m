%%% PROGRAM TO OBTAIN FFT PLOT OF RAW DATA OBTAINED USING DCA1000 AND AWR1243 
%% Only adc_data.bin is required. 
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
fid = fopen('adc_data_Raw_Raw_0.bin','r');
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


 epsilon = 0.11;
 epsilon_dash = 0.06;
 delta = 2*10^-9;
 k = 20;
 n = 256;
 
 % Window array shape must be 1x256
    % Dolph Chebyshev Window convolved with box car
    H = sigwin.chebwin(256,200);
    win = generate(H);
    WIN = fft(win);
    box_car = fft(ones(round(2*epsilon_dash*n),1),256);
    abs_win = abs(WIN);
    final_win = box_car .* (abs(WIN)/max(abs_win));
    g= transpose(ifft(final_win,256));
    g = 4.38*(g/max(g));
    
    g(1:128) = flip(g(129:256));
    G = fft(g,256);
    
    % Hanning Window
    % g= transpose(hanning(256));
    
    % Bartlett Window
    % g= transpose(bartlett(256));

    % Blackman Window
    % g= transpose(blackman(256));

    % Hamming Window
    % g= transpose(hamming(256));

    % Taylor Window
    % g= transpose(taylorwin(256));

    % Kaiser Window
    % g= transpose(kaiser(256,2.5));

    % disp(size(g));
    
    
    % plot(0:255,20*log10(1+abs(fft(g))));
    % title('FFT of window');
    % ylabel('FFT [dB]');
    % figure;

    % plot(0:255,abs(fft(g)));
    % title('FFT of window');
    % ylabel('FFT');
    % figure;
    
    % plot(0:255,g);
    % title('Window');
    % figure;

    % plot(0:255,log10(1+g));
    % title('Window');
    % ylabel('dB');
    % figure;


    % Loop for finding the b value
    % for c = 0.2:0.01:1
    %     B = round((c*k)/epsilon);
    %     if  mod(n,B) == 0 
    %         break
    %     end
    % end
    
    B =32; % or change it to 8

    d = round(0.1/epsilon);
    
    % d = 1;

    disp("B");
    disp(B);
    disp("d");
    disp(d);
    runtime_sfft = 0;
    runtime_inbuilt = 0;

%  for  idx=1:Dx:plotEnd 
    % for  idx=1:Dx:plotEnd/3 

for  idx=1:Dx:Dx 
    clf;


    R1=real_1(idx : idx+Dx-1);
    I1 = imag_1(idx : idx+Dx-1);
    
    tic
    R11 = fftshift(fft_recur(hanning(length(R1)).*R1));
    I11 = fftshift(fft_recur(hanning(length(I1)).*I1));
    ABS11 = fftshift(fft_recur(hanning(length(R1+1i*I1)).*(R1+1i*I1)));
  
    runtime_inbuilt = runtime_inbuilt + toc;
       
    x = R1 + 1i*I1;
    n = length(x);  
    

    [I,h_sigma,o_sigma,Z,tau] = inner_location_loop(x,n,k,g,d,B,delta);
    
    tic
    sfft_real = fftshift(outer_loop(R1,n,k,g,d,B,delta,G));
    sfft_imag = fftshift(outer_loop(I1,n,k,g,d,B,delta,G));
    sfft = fftshift(outer_loop(R1+1i*I1,n,k,g,d,B,delta,G));
    % sfft_real = sfft_real + flip(sfft_real)
    % sfft_imag = sfft_imag + flip(sfft_imag)
    % sfft = sfft + flip(sfft)
    runtime_sfft = runtime_sfft + toc;
    
%%%%%%%%%% TO RUN BIGBAND uncomment the following lines %%%%%%%%%%%%%%%
% 	R1_dash = R1(1:255)
% 	I1_dash = I1(1:255)
% 	x_dash = x(1:255)
% coprimes = [3,5];

% 	   sfft_real  = bigband(R1_dash,coprimes,255)
% 	   sfft_imag  = bigband(I1_dash,coprimes,255)
% 	   sfft  = bigband(x_dash,coprimes,255)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
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
%     axis([-25 25 -140 0]);
   title('SFFT Amplitude(per chirp-real only)','FontSize',20);
   xlabel('Distance[m]','FontSize',18);
   ylabel('FFT output [dBFS]','FontSize',18);
   set(gca,'FontSize',18,'FontWeight','bold');

   
  subplot(2,2,3);
    grid on ;
    hold all;
    stem(distance,dBFS_imag_s,'linewidth',2);
    % axis([-25 25 -140 0]);
%     axis([-25 25 -140 0]);
    title('SFFT Amplitude(per chirp-imag only)','FontSize',20);
    xlabel('Distance[m]','FontSize',18);
    ylabel('FFT output [dBFS]','FontSize',18);
    set(gca,'FontSize',18,'FontWeight','bold');

   subplot(2,2,4);
    grid on ;
    hold all;
    stem(distance,dBFS_cmplx_s,'linewidth',1);
    % axis([-25 25 -140 0]);
%     axis([-25 25 -140 0]);
    title('SFFT Amplitude(per chirp)','FontSize',20);
    xlabel('Distance[m]','FontSize',18);
    ylabel('FFT output [dBFS]','FontSize',18);
    set(gca,'FontSize',18,'FontWeight','bold');
   
    % drawnow;
   
    
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
%     axis([-25 25 -140 0]);
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
%     axis([-25 25 -140 0]);
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
%     axis([-25 25 -140 0]);
     title('1D FFT Amplitude profile(per chirp)','FontSize',20);
     xlabel('Distance[m]','FontSize',18);
     ylabel('FFT output [dBFS]','FontSize',18);
     set(gca,'FontSize',18,'FontWeight','bold');
     legend('SFFT','Normal FFT','FontSize',7)
     drawnow;

% indices = (dBFS_real_s ~= -Inf);
% fprintf("  Expected Sparsity : %d \n ",k);
% fprintf("  Calculated Sparsity %d  \n \n",sum(indices));
% err = sum((dBFS_real(indices)-dBFS_real_s(indices)).^2)/length(indices);
% r_mean = sum(dBFS_real(indices)./dBFS_real_s(indices))/length(indices);
% fprintf(" MEAN SQUARE ERROR Loop index %d in real values is %d \n",int64((idx+Dx-1))/Dx,err);
% fprintf(" MEAN of Value ratio Loop index %d in real values is %d \n \n",int64((idx+Dx-1))/Dx,r_mean);
% indices = (dBFS_imag_s ~= -Inf);
% fprintf("  Expected Sparsity : %d \n ",k);
% fprintf("  Calculated Sparsity %d  \n \n",sum(indices));
% err = sum((dBFS_imag(indices)-dBFS_imag_s(indices)).^2)/length(indices);
% r_mean = sum(dBFS_imag(indices)./dBFS_imag_s(indices))/length(indices);
% fprintf(" MEAN SQUARE ERROR Loop index %d in imag values is %d \n",int64((idx+Dx-1)/Dx),err);
% fprintf(" MEAN of Value ratio Loop index %d in imag values is %d \n \n",int64((idx+Dx-1))/Dx,r_mean);
% indices = (dBFS_cmplx_s ~= -Inf);
% fprintf("  Expected Sparsity : %d \n ",k);
% fprintf("  Calculated Sparsity %d  \n \n",sum(indices));
% err = sum((dBFS_cmplx(indices)-dBFS_cmplx_s(indices)).^2)/length(indices);
% r_mean = sum(dBFS_cmplx(indices)./dBFS_cmplx_s(indices))/length(indices);
% fprintf(" MEAN SQUARE ERROR Loop index %d in cmplx values is %d \n ",int64((idx+Dx-1))/Dx,err);
% fprintf(" MEAN of Value ratio Loop index %d in cmplx values is %d \n \n",int64((idx+Dx-1))/Dx,r_mean);

% indices = (abs(sfft_real) ~= 0);
% err = sum((R11(indices)-sfft_real(indices)).^2)/length(indices);
% r_mean = sum(R11(indices)./sfft_real(indices))/length(indices);
% fprintf(" MEAN SQUARE ERROR Loop index %d in real values is %d \n",int64((idx+Dx-1))/Dx,err);
% fprintf(" MEAN of Value ratio Loop index %d in real values is %d \n \n",int64((idx+Dx-1))/Dx,r_mean);
% indices = (abs(sfft_imag) ~= 0);
% err = sum((I11(indices)-sfft_imag(indices)).^2)/length(indices);
% r_mean = sum(I11(indices)./sfft_imag(indices))/length(indices);
% fprintf(" MEAN SQUARE ERROR Loop index %d in imag values is %d \n",int64((idx+Dx-1)/Dx),err);
% fprintf(" MEAN of Value ratio Loop index %d in imag values is %d \n \n",int64((idx+Dx-1))/Dx,r_mean);
% indices = (abs(sfft) ~= 0);
% err = sum((ABS11(indices)-sfft(indices)).^2)/length(indices);
% r_mean = sum(ABS11(indices)./sfft(indices))/length(indices);
% fprintf(" MEAN SQUARE ERROR Loop index %d in cmplx values is %d \n ",int64((idx+Dx-1))/Dx,err);
% fprintf(" MEAN of Value ratio Loop index %d in cmplx values is %d \n \n",int64((idx+Dx-1))/Dx,r_mean);


    %  pause(2);


    

end
fprintf("  Inbuilt Runtime = %d \n",runtime_inbuilt);
fprintf("  Sfft Runtime = %d \n",runtime_sfft);

    % 
% % plotting for all frames (Channel 2)
% 
real_2 = adcData(:,2);         
imag_2 = adcData(:,6);
% 
% figure(2);
% hold on;


%  for  idx=1:Dx:plotEnd 
% for  idx=1:Dx:plotEnd/3 
% % for  idx=1:Dx:Dx 
%     clf;    
%     R2=real_2(idx : idx+Dx-1);
%     I2 = imag_2(idx : idx+Dx-1);




%    R22 = fftshift(fft(hanning(length(R2)).*R2));
%    I22 = fftshift(fft(hanning(length(I2)).*I2));
%   ABS22 = fftshift(fft(hanning(length(R2+1i*I2)).*(R2+1i*I2)));
 

%     sfft_real = fftshift(outer_loop(R2,n,k,g,d,B,delta,G));
%     sfft_imag = fftshift(outer_loop(I2,n,k,g,d,B,delta,G));
%     sfft = fftshift(outer_loop(R2+1i*I2,n,k,g,d,B,delta,G));
    
    
%     dB_cmplx = 20*log10(abs(sfft));  
%     dB_real = 20*log10(abs(sfft_real));
%     dB_imag = 20*log10(abs(sfft_imag));
%     dBFS_real_s = dB_real - 20*log10(256)-20*log10(2^15)+20*log10(2^(0.5))-2.0;
%     dBFS_imag_s = dB_imag - 20*log10(256)-20*log10(2^15)+20*log10(2^(0.5))-2.0;
%     dBFS_cmplx_s = dB_cmplx - 20*log10(256)-20*log10(2^15)+20*log10(2^(0.5))-2.0;
    

%  subplot(2,2,1)  
%  hold all;
%     plot( distance,R2,'b','linewidth',1);
%     plot( distance,I2,'r','linewidth',1);
%     title('IQ time domain signal','FontSize',20);
%     axis([-25 25 -500 500]);
%     xlabel('Samples','FontSize',18);
%     ylabel('ADC value','FontSize',18);
%     grid on ;
%     set(gca,'FontSize',18,'FontWeight','bold');
     

%   subplot(2,2,2);  
%    grid on ;
%    hold all;
%    stem(distance,dBFS_real_s,'linewidth',2);
%    axis([-25 25 -140 0]);
%    title('1D FFT Amplitude profile(per chirp-real only)','FontSize',20);
%    xlabel('Distance[m]','FontSize',18);
%    ylabel('FFT output [dBFS]','FontSize',18);
%    set(gca,'FontSize',18,'FontWeight','bold');

   
%   subplot(2,2,3);
%     grid on ;
%     hold all;
%     stem(distance,dBFS_imag_s,'linewidth',2);
%     axis([-25 25 -140 0]);
%     title('1D FFT Amplitude profile(per chirp-imag only)','FontSize',20);
%     xlabel('Distance[m]','FontSize',18);
%     ylabel('FFT output [dBFS]','FontSize',18);
%     set(gca,'FontSize',18,'FontWeight','bold');

%    subplot(2,2,4);
%     grid on ;
%     hold all;
%     stem(distance,dBFS_cmplx_s,'linewidth',1);
%     axis([-25 25 -140 0]);
%     title('1D FFT Amplitude profile(per chirp)','FontSize',20);
%     xlabel('Distance[m]','FontSize',18);
%     ylabel('FFT output [dBFS]','FontSize',18);
%     set(gca,'FontSize',18,'FontWeight','bold');
   
%   dB_cmplx = 20*log10(abs(ABS22));
%   dB_real = 20*log10(abs(R22));
%   dB_imag = 20*log10(abs(I22));
%   dBFS_real = dB_real - 20*log10(256)-20*log10(2^15)+20*log10(2^(0.5));
%   dBFS_imag = dB_imag - 20*log10(256)-20*log10(2^15)+20*log10(2^(0.5));
%   dBFS_cmplx= dB_cmplx - 20*log10(256)-20*log10(2^15)+20*log10(2^(0.5));
   
    
%     % figure;
%     subplot(2,2,1)
%     hold all;
%     plot( distance,R2,'b','linewidth',1);
%     plot( distance,I2,'r','linewidth',1);
%     title('IQ time domain signal','FontSize',20);
%     axis([-25 25 -500 500]);
%     xlabel('Samples','FontSize',18);
%     ylabel('ADC value','FontSize',18);
%     grid on ;
%     set(gca,'FontSize',18,'FontWeight','bold');


%     subplot(2,2,2);  
%     grid on ;
%     hold all;
%     plot(distance,dBFS_real,'k','linewidth',2);
%     axis([-25 25 -140 0]);
%     %     axis([-25 25 -140 0]);
%     title('1D FFT Amplitude profile(per chirp-real only)','FontSize',20);
%     xlabel('Distance[m]','FontSize',18);
%     ylabel('FFT output [dBFS]','FontSize',18);
%     set(gca,'FontSize',18,'FontWeight','bold');
%     legend('SFFT','Normal FFT','FontSize',7)


%     subplot(2,2,3);
%     grid on ;
%     hold all;
%     plot(distance,dBFS_imag,'g','linewidth',2);
%     axis([-25 25 -140 0]);
%     %     axis([-25 25 -140 0]);
%     title('1D FFT Amplitude profile(per chirp-imag only)','FontSize',20);
%     xlabel('Distance[m]','FontSize',18);
%     ylabel('FFT output [dBFS]','FontSize',18);
%     set(gca,'FontSize',18,'FontWeight','bold');
%     legend('SFFT','Normal FFT','FontSize',7)

%     subplot(2,2,4);
%     grid on ;
%     hold all;
%     plot(distance,dBFS_cmplx,'r','linewidth',1);
%     axis([-25 25 -140 0]);
%     %     axis([-25 25 -140 0]);
%     title('1D FFT Amplitude profile(per chirp)','FontSize',20);
%     xlabel('Distance[m]','FontSize',18);
%     ylabel('FFT output [dBFS]','FontSize',18);
%     set(gca,'FontSize',18,'FontWeight','bold');
%     drawnow;
%     legend('SFFT','Normal FFT','FontSize',7)

%    drawnow;

%     indices = (dBFS_real_s ~= -Inf);
%     % disp(indices)
%     err = sum((dBFS_real(indices)-dBFS_real_s(indices)).^2)/length(indices);
%     r_mean = sum(dBFS_real(indices)./dBFS_real_s(indices))/length(indices);
%     % fprintf(" MEAN SQUARE ERROR Loop index %d in real values is %d \n",int64((idx+Dx-1))/Dx,err);
%     % fprintf(" MEAN of Value ratio index %d in real values is %d \n \n",int64((idx+Dx-1))/Dx,r_mean);
%     indices = (dBFS_imag_s ~= -Inf);
%     err = sum((dBFS_imag(indices)-dBFS_imag_s(indices)).^2)/length(indices);
%     r_mean = sum(dBFS_imag(indices)./dBFS_imag_s(indices))/length(indices);
%     % fprintf(" MEAN SQUARE ERROR Loop index %d in imag values is %d \n",int64((idx+Dx-1)/Dx),err);
%     % fprintf(" MEAN of Value ratio index %d in imag values is %d \n \n",int64((idx+Dx-1))/Dx,r_mean);
%     indices = (dBFS_cmplx_s ~= -Inf);
%     err = sum((dBFS_cmplx(indices)-dBFS_cmplx_s(indices)).^2)/length(indices);
%     r_mean = sum(dBFS_cmplx(indices)./dBFS_cmplx_s(indices))/length(indices);
%     % fprintf(" MEAN SQUARE ERROR Loop index %d in cmplx values is %d \n \n",int64((idx+Dx-1))/Dx,err);
%     % fprintf(" MEAN of Value ratio index %d in cmplx values is %d \n \n",int64((idx+Dx-1))/Dx,r_mean);


%     %  pause(2);

% end
% 
% 
% % plotting for all frames (Channel 3)
% 
 real_3 = adcData(:,3);         
imag_3 = adcData(:,7);
% 
% figure(3);                    
% hold on ;
% for  n=1:Dx:plotEnd 
%     clf;
%    R3=real_3(n : n+Dx-1);
%    I3 = imag_3(n : n+Dx-1);
%    R33 =fftshift(fft(hanning(length(R3)).*R3));
%    I33 =fftshift(fft(hanning(length(I3)).*I3));
%   ABS33 = fftshift(fft(hanning(length(R3+1i*I3)).*(R3+1i*I3)));
% 
%   dB_cmplx = 20*log10(abs(ABS33));
%   dB_real = 20*log10(abs(R33));
%   dB_imag = 20*log10(abs(I33));
%   dBFS_real = dB_real - 20*log10(256)-20*log10(2^15)+20*log10(2^(0.5));
%   dBFS_imag = dB_imag - 20*log10(256)-20*log10(2^15)+20*log10(2^(0.5));
%   dBFS_cmplx= dB_cmplx - 20*log10(256)-20*log10(2^15)+20*log10(2^(0.5));
% 
%  subplot(2,2,1)  
%  hold all;
%     plot( distance,R3,'b','linewidth',1);
%     plot( distance,I3,'r','linewidth',1);
%     title('IQ time domain signal','FontSize',20);
%     axis([-25 25 -500 500]);
%     xlabel('Samples','FontSize',18);
%     ylabel('ADC value','FontSize',18);
%     grid on ;
%     set(gca,'FontSize',18,'FontWeight','bold');
%      
%  %  subplot(2,3,2)
% %     plot( R3,'b','linewidth',2);
% %     title('IQ time domain signal for real data','FontSize',20);
% %     xlabel('time(*0.1 us)','FontSize',18);
% %     ylabel('adc value','FontSize',18);
% %    axis([-25 25 -500 500]);
% %      grid on ;
% 
% %   subplot(2,3,3)
% %    plot( I3,'r','linewidth',2); 
% %    title('IQ time domain signal for imaginary data','FontSize',20);
%    xlabel('time(*0.1 us)','FontSize',18);
%    ylabel('adc value','FontSize',18);
%    axis([-25 25 -500 500]);
%    grid on ;
%    set(gca,'FontSize',18,'FontWeight','bold')

%   subplot(2,2,2);  
%    grid on ;
%    hold all;
%    plot(distance,dBFS_real,'b','linewidth',2);
%    axis([-25 25 -140 0]);
%    title('1D FFT Amplitude profile(per chirp-real only)','FontSize',20);
%    xlabel('Distance[m]','FontSize',18);
%    ylabel('FFT output [dBFS]','FontSize',18);
%    set(gca,'FontSize',18,'FontWeight','bold');
% 
%    
%   subplot(2,2,3);
%     grid on ;
%     hold all;
%     plot(distance,dBFS_imag,'g','linewidth',2);
%     axis([-25 25 -140 0]);
%     title('1D FFT Amplitude profile(per chirp-imag only)','FontSize',20);
%     xlabel('Distance[m]','FontSize',18);
%     ylabel('FFT output [dBFS]','FontSize',18);
%     set(gca,'FontSize',18,'FontWeight','bold');
% 
%    subplot(2,2,4);
%     grid on ;
%     hold all;
%     plot(distance,dBFS_cmplx,'r','linewidth',1);
%     axis([-25 25 -140 0]);
%     title('1D FFT Amplitude profile(per chirp)','FontSize',20);
%     xlabel('Distance[m]','FontSize',18);
%     ylabel('FFT output [dBFS]','FontSize',18);
%     set(gca,'FontSize',18,'FontWeight','bold');
%    drawnow;
% end
% 
% % plotting for all frames (Channel 4)
% 
 real_4 = adcData(:,4);         
imag_4 = adcData(:,8);
% 
% figure(4);                    
% hold on ;
% for  n=1:Dx:plotEnd 
%     clf;
%    R4=real_4(n : n+Dx-1);
%    I4 = imag_4(n : n+Dx-1);
%    R44 = fftshift(fft(hanning(length(R4)).*R4));
%    I44 = fftshift(fft(hanning(length(I4)).*I4));
%   ABS44 = fftshift(fft(hanning(length(R4+1i*I4)).*(R4+1i*I4)));
% 
%   dB_cmplx = 20*log10(abs(ABS44));
%   dB_real = 20*log10(abs(R44));
%   dB_imag = 20*log10(abs(I44));
%   dBFS_real = dB_real - 20*log10(256)-20*log10(2^15)+20*log10(2^(0.5));
%   dBFS_imag = dB_imag - 20*log10(256)-20*log10(2^15)+20*log10(2^(0.5));
%   dBFS_cmplx= dB_cmplx - 20*log10(256)-20*log10(2^15)+20*log10(2^(0.5));
% 
%  subplot(2,2,1)  
%  hold all;
%     plot( distance,R4,'b','linewidth',1);
%     plot( distance,I4,'r','linewidth',1);
%     title('IQ time domain signal','FontSize',20);
%     axis([-25 25 -500 500]);
%     xlabel('Samples','FontSize',18);
%     ylabel('ADC value','FontSize',18);
%     grid on ;
%     set(gca,'FontSize',18,'FontWeight','bold');
%      
%  %  subplot(2,3,2)
% %     plot( R4,'b','linewidth',2);
% %     title('IQ time domain signal for real data','FontSize',20);
% %     xlabel('time(*0.1 us)','FontSize',18);
% %     ylabel('adc value','FontSize',18);
% %    axis([-25 25 -500 500]);
% %      grid on ;
% 
% %   subplot(2,3,3)
% %    plot( I4,'r','linewidth',2); 
% %    title('IQ time domain signal for imaginary data','FontSize',20);
% %    xlabel('time(*0.1 us)','FontSize',18);
% %    ylabel('adc value','FontSize',18);
% %    axis([-25 25 -500 500]);
% %    grid on ;
% %    set(gca,'FontSize',18,'FontWeight','bold')
% 
%   subplot(2,2,2);  
%    grid on ;
%    hold all;
%    plot(distance,dBFS_real,'b','linewidth',2);
%    axis([-25 25 -140 0]);
%    title('1D FFT Amplitude profile(per chirp-real only)','FontSize',20);
%    xlabel('Distance[m]','FontSize',18);
%    ylabel('FFT output [dBFS]','FontSize',18);
%    set(gca,'FontSize',18,'FontWeight','bold');
% 
%    
%   subplot(2,2,3);
%     grid on ;
%     hold all;
%     plot(distance,dBFS_imag,'g','linewidth',2);
%     axis([-25 25 -140 0]);
%     title('1D FFT Amplitude profile(per chirp-imag only)','FontSize',20);
%     xlabel('Distance[m]','FontSize',18);
%     ylabel('FFT output [dBFS]','FontSize',18);
%     set(gca,'FontSize',18,'FontWeight','bold');
% 
%    subplot(2,2,4);
%     grid on ;
%     hold all;
%     plot(distance,dBFS_cmplx,'r','linewidth',1);
%     axis([-25 25 -140 0]);
%     title('1D FFT Amplitude profile(per chirp)','FontSize',20);
%     xlabel('Distance[m]','FontSize',18);
%     ylabel('FFT output [dBFS]','FontSize',18);
%     set(gca,'FontSize',18,'FontWeight','bold');
%    drawnow;
% end
% 
% FFT of all the channels for comparison

% figure(5);
% hold on;
% for  n=1:Dx:plotEnd 
%     clf;
%      R1=real_1(n : n+Dx-1);
%      I1 = imag_1(n : n+Dx-1);
%      ABS11 = fftshift(fft(hanning(length(R1+1i*I1)).*(R1+1i*I1)));
%      dB_cmplx_1 = 20*log10(abs(ABS11));  
%      dBFS_cmplx_1= dB_cmplx_1 - 20*log10(256)-20*log10(2^15)+20*log10(2^(0.5))-2.0;

%      R2=real_2(n : n+Dx-1);
%      I2 = imag_2(n : n+Dx-1);
%      ABS22 = fftshift(fft(hanning(length(R2+1i*I2)).*(R2+1i*I2)));
%      dB_cmplx_2 = 20*log10(abs(ABS22));  
%      dBFS_cmplx_2= dB_cmplx_2 - 20*log10(256)-20*log10(2^15)+20*log10(2^(0.5))-2.0;

  
%      R3=real_3(n : n+Dx-1);
%      I3 = imag_3(n : n+Dx-1);
%      ABS33 = fftshift(fft(hanning(length(R3+1i*I3)).*(R3+1i*I3)));
%      dB_cmplx_3 = 20*log10(abs(ABS33));  
%      dBFS_cmplx_3= dB_cmplx_3 - 20*log10(256)-20*log10(2^15)+20*log10(2^(0.5))-2.0;
   
%      R4=real_4(n : n+Dx-1);
%      I4 = imag_4(n : n+Dx-1);
%      ABS44 = fftshift(fft(hanning(length(R4+1i*I4)).*(R4+1i*I4)));
%      dB_cmplx_4 = 20*log10(abs(ABS44));  
%      dBFS_cmplx_4= dB_cmplx_4 - 20*log10(256)-20*log10(2^15)+20*log10(2^(0.5))-2.0;
   
%  subplot(2,2,1);
%   grid on;
%   hold all;
%   axis([-25 25 -140 0]);
%   title('FFT of Channel-1 signal with mirror (dB)');
%   xlabel('Distance[m]','FontSize',18);
%   ylabel('FFT output [dBFS]','FontSize',18);
%   set(gca,'FontSize',18,'FontWeight','bold');
%   plot(distance,dBFS_cmplx_1);

% subplot(2,2,2);
%   grid on;
%   hold all;
%   axis([-25 25 -140 0]);
%   title('FFT of Channel-2 signal with mirror (dB)');
%   xlabel('Distance[m]','FontSize',18);
%   ylabel('FFT output [dBFS]','FontSize',18);
%   set(gca,'FontSize',18,'FontWeight','bold');
%   plot(distance,dBFS_cmplx_2);

% subplot(2,2,3);
%   grid on;
%   hold all;
%   axis([-25 25 -140 0]);
%   title('FFT of Channel-3 signal with mirror (dB)');
%   xlabel('Distance[m]','FontSize',18);
%   ylabel('FFT output [dBFS]','FontSize',18);
%   set(gca,'FontSize',18,'FontWeight','bold');
%   plot(distance,dBFS_cmplx_3);

% subplot(2,2,4);
%   grid on;
%   hold all;
%   axis([-25 25 -140 0]);
%   title('FFT of Channel-4 signal with mirror (dB)');
%   xlabel('Distance[m]','FontSize',18);
%   ylabel('FFT output [dBFS]','FontSize',18);
%   set(gca,'FontSize',18,'FontWeight','bold');
%   plot(distance,dBFS_cmplx_4);
% drawnow;
% end


