%
% Computing the SER for M-QAM, M-PSK, M-PAM modulations
% in a Rayleigh fading channel with Alamouti Space Time Block Coding(STBC)
% transmit diversity scheme
%
% Coded by Amir otd :)[https://github.com/amirotd]
%

clc; clear; close all;

alamouti(16, 'QAM');
alamouti(8, 'PSK');
alamouti(4, 'PAM');

function alamouti(M, mod_type)
    N = 1e5;
    snr_db = 0:1:20;
    
    tx_data = randi([0 M-1],1,2*N);
    temp1 = reshape(tx_data,N,2);
    
    switch mod_type
        % Modulation
        case{'PSK'}
            moded = pskmod(temp1,M);
        case{'QAM'}
            moded = qammod(temp1,M);
        case{'PAM'}
            moded = pammod(temp1,M);
    end

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
        
        switch mod_type
            % Demodulation
            case{'PSK'}
                rx_data(:,1) = pskdemod(z1(:,1),M);
                rx_data(:,2) = pskdemod(z1(:,2),M);
            case{'QAM'}
                rx_data(:,1) = qamdemod(z1(:,1),M);
                rx_data(:,2) = qamdemod(z1(:,2),M);
            case{'PAM'}
                rx_data(:,1) = pamdemod(z1(:,1),M);
                rx_data(:,2) = pamdemod(z1(:,2),M);
        end

        % Get the Symbol Error Rate   
        [~, alam_ser(i)] = symerr(rx_data, temp1);
        
        %% Without Alamouti
        y2 = awgn(x2,snr_db(i),'measured'); % Added White Gaussian noise
        z2 = y2./h;
    
        switch mod_type
            % Demodulation
            case{'PSK'}
                rx_data2 = pskdemod(z2,M);
            case{'QAM'}
                rx_data2 = qamdemod(z2,M);
            case{'PAM'}
                rx_data2 = pamdemod(z2,M);
        end
        
        % Get the Symbol Error Rate
        [~, no_alam_ser(i)] = symerr(rx_data2, temp1);
    end
    
    % plot SER vs SNR
    figure
    semilogy(snr_db,alam_ser,'r*-');
    hold on
    semilogy(snr_db,no_alam_ser,'bo-');
    xlabel('SNR (dB)');
    ylabel('Symbol Error Rate (SER)');
    title(['Alamouti SER for ', num2str(M), mod_type, ' Modulation'])
    legend('With Alamouti (2x1 channel)','Without Alamouti (1x1 channel)')
    grid on
end
