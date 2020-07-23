function [Capacity] = MIMO_Capacity(mode, Nr, Nt)
% Description:
% ------------- 
% Calculate the ergodic capacity of a MIMO system using optimal power
% allocation of water filling, given that channel state information is 
% available at both the transmitter and receiver.
% Syntax: 
% ------------
% function MIMO_Capacity(fid, variable, value, units)
%
% Inputs:
% ---------
%    mode - 0: assume there is no correlation at Tx or Rx; 
%              the efficiency of all antennas is 100%
%           1: assume correlation exists at Tx only;
%              only transmit antennas are lossy
%           2: assume correlation exists at both the Tx and Rx;
%              all antennas are lossy
%
%    Nr - Number of receive antennas
%    Nt - Number of transmit antennas
%    
%
% Outputs:
% ---------
%    Capacity - Ergodic capacity of the MIMO system
%
% Example: 
% ---------
%    Line 1 of example
%    Line 2 of example
%    Line 3 of example
%
% Other m-files required: none
% Subfunctions: function [C_for_one_iter] = Water_Filling(Pt, Noise_Power, Nr, Nt, R_Rx, R_Tx)
% MAT-files required: none
%
% See also: OTHER_FUNCTION_NAME1,  OTHER_FUNCTION_NAME2
% ----------------------------------------------------------------------------
% Author: Dr. LIU Yanting
% State Key Laboratory of Terahertz and Millimeter Waves, City University of Hong Kong
% Email address: yantliu3-c@my.cityu.edu.hk
% Jul.2020
% ----------------------------------------------------------------------------

%------------- BEGIN CODE --------------

if nargin == 1
    Nr = 3;
    Nt = 3;
end

N_iter = 10000; % Number of iterations
Noise_Power = 1; % Noise power sigma_n^2
Pt = 100; % Transmit power

% Empty matrix for storing capacity in each iteration
C = zeros(N_iter, 1);

switch mode
    case 0
        for iter = 1:N_iter
                C(iter) = Water_Filling(Pt, Noise_Power, Nr, Nt);
        end
        Capacity = mean(C); % Mean of all capacities
    
    case {1,2}
        ECC = readtable('ECC.xlsx');
        % Empty matrix for storing mean capacity over all iterations
        % at each frequency point
        Capacity = zeros(length(ECC.Freq),1);
        
        % Load ECC and efficiency data
        for freq_iter = 1:length(ECC.Freq)
            Corr_Tx = eye(3);
            Corr_Tx(1,2) = ECC.P12(freq_iter);
            Corr_Tx(2,1) = ECC.P12(freq_iter);
            Corr_Tx(1,3) = ECC.P13(freq_iter);
            Corr_Tx(3,1) = ECC.P13(freq_iter);
            Corr_Tx(2,3) = ECC.P23(freq_iter);
            Corr_Tx(3,2) = ECC.P23(freq_iter);
            % Cholesky factorization
            R_Tx = chol(Corr_Tx);
            
            if mode == 2
                R_Rx = R_Tx;
                Pt_real = Pt * (ECC.Eff(freq_iter))^2;
            else
                R_Rx = eye(3);
                Pt_real = Pt * ECC.Eff(freq_iter);
            end
            

            for iter = 1:N_iter
                C(iter) = Water_Filling(Pt_real, Noise_Power, Nr, Nt, R_Rx, R_Tx);
            end

            Capacity(freq_iter) = mean(C); % Mean of all capacities
        end

end
end


function [C_for_one_iter] = Water_Filling(Pt, Noise_Power, Nr, Nt, R_Rx, R_Tx)
   
   % Channel matrix given that all paths are i.i.d. Gaussian distributions
   % with a zero mean and a variance of 1/2 i.e. hij ~ N(0, 1/2) 
    H = sqrt(1/2)*(randn(Nr,Nt)+1i*randn(Nr,Nt));
   
   % If number of function input arguments is 6, it should be mode 1 or
   % mode 2 and channel matrix should be the product of transmit matrix,
   % Rayleigh distribution matrix and receive matrix.
    if nargin == 6
        H = R_Rx * H * R_Tx;
    end
    
    % Square of singular values which represent the channel gains of
    % spatial multiplexing
    sv_square = abs(svd(H)).^2; 
    
    % Screen out trivial values
    sv_square = sv_square(sv_square>0);
    
    % Number of nontrivial singular values
    N_sv = length(sv_square);
    % Ratio of noise power to channel gains
    r1 = Noise_Power./sv_square;

    % Lagrange multiplier
    lambda = N_sv/((sum(r1)+Pt)*log(2));
    
    % Minimum allocated power using water filling scheme
    threshold = 1/(lambda*log(2)) - r1(end);

    % Check if the minimum allocated power is <= 0. If so, screen out
    % this "bad" channel and allocate transmit power to the remaining 
    % parallel channels. Repeat this process until all allocated power 
    % is positive.  
     while 1
        if threshold <= 0
            sv_square = sv_square(1:end-1);
            N_sv = N_sv - 1;
            r1 = Noise_Power./sv_square;
            lambda = N_sv/((sum(r1)+Pt)*log(2));
            threshold = 1/(lambda*log(2)) - r1(end);
        else
            break
        end
     end

    % Final optimal power allocation scheme  
    Power_allocation = 1/(lambda*log(2)) - r1;
    
    % Total normalized Shannon capacity (bit/s/Hz)
    C_for_one_iter = sum(log2(1 + Power_allocation./r1)); 
end


%------------- END OF CODE --------------