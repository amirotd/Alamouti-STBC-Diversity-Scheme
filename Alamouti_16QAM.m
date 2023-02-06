%
% Computing the SER for 16QAM modulation in a Rayleigh fading channel
% with Alamouti Space Time Block Coding(STBC)
% transmit diversity scheme
%
% Coded by Amir otd :)[https://github.com/amirotd]
%

clc; clear; close all;

N = 1e5;
snr_db = 0:1:20;

% Generate random 16QAM symbols
alpha16qam = [-3 -1 1 3];
tx_data = randsrc(1,2*N,alpha16qam) + 1i*randsrc(1,2*N,alpha16qam);
moded = reshape(tx_data,N,2);

% Rayleigh fading channel
h  = 1/sqrt(2) * (randn(1,2*N) + 1i*randn(1,2*N));
h = reshape(h,N,2);

% TX-Data over a Rayleigh fading channel in alamouti scheme
x1(:,1) =  h(:,1).* 1/sqrt(2).*moded(:,1) + h(:,2).* 1/sqrt(2).*moded(:,2);
x1(:,2) = -h(:,1).*conj(1/sqrt(2).*moded(:,2)) + h(:,2).*conj(1/sqrt(2).*moded(:,1));

% TX-Data over a Rayleigh fading channel(without alamouti)
x2 = h.* ((1/sqrt(2)).*moded);

for i = 1 : length(snr_db)
    
    %% With Alamouti
    y1 = awgn(x1,snr_db(i),'measured'); % Added White Gaussian noise
    
    z1(:,1) = sqrt(2)*(conj(h(:,1)) .* y1(:,1) + h(:,2) .* conj(y1(:,2)))./(abs(h(:,1)).^2 + abs(h(:,2)).^2);
    z1(:,2) = sqrt(2)*(conj(h(:,2)) .* y1(:,1) - h(:,1) .* conj(y1(:,2)))./(abs(h(:,1)).^2 + abs(h(:,2)).^2);
    
    % Demodulation process
    rx_data_re((real(z1)< -2)) = -3;
    rx_data_re((real(z1) > 2)) =  3;
    rx_data_re((real(z1)>-2 & real(z1)<=0)) = -1;
    rx_data_re((real(z1)>0 & real(z1)<=2)) = 1;
    
    rx_data_im((imag(z1)< -2)) = -3;
    rx_data_im((imag(z1) > 2)) =  3;
    rx_data_im((imag(z1)>-2 & imag(z1)<=0)) = -1;
    rx_data_im((imag(z1)>0 & imag(z1)<=2)) = 1;
    rx_data = rx_data_re + 1i*rx_data_im;

    % Get the Symbol Error Rate
    [~, alam_ser(i)] = symerr(rx_data, tx_data);
    
    %% Without Alamouti
    y2 = awgn(x2,snr_db(i),'measured'); % Added White Gaussian noise
    z2 = y2./h;
    
    % Demodulation process
    rx_data_re2((real(z2)< -2)) = -3;
    rx_data_re2((real(z2) > 2)) =  3;
    rx_data_re2((real(z2)>-2 & real(z2)<=0)) = -1;
    rx_data_re2((real(z2)>0 & real(z2)<=2)) = 1;
    
    rx_data_im2((imag(z2)< -2)) = -3;
    rx_data_im2((imag(z2) > 2)) =  3;
    rx_data_im2((imag(z2)>-2 & imag(z2)<=0)) = -1;
    rx_data_im2((imag(z2)>0 & imag(z2)<=2)) = 1;
    rx_data2 = rx_data_re2 + 1i*rx_data_im2;
    
    % Get the Symbol Error Rate
    [~, no_alam_ser(i)] = symerr(rx_data2, tx_data);

end

% plot SER vs SNR
semilogy(snr_db,alam_ser,'r*-');
hold on
semilogy(snr_db,no_alam_ser,'bo-');
xlabel('SNR (dB)');
ylabel('Symbol Error Rate (SER)');
title('Alamouti SER for 16QAM Modulation')
legend('With Alamouti (2x1 channel)','Without Alamouti (1x1 channel)')
grid on
