% Serially Concatenated Convolutional Coded (SCCC) QPSK over AWGN - logMAP

close all
clear all
clc
%---------------- SIMULATION PARAMETERS ------------------------------------
SNR_dB = 1.5; % SNR per bit in dB (in logarithmic scale)
sim_runs = 1*(10^1); % simulation runs
frame_size = 1024; % frame size
num_bit = 0.5*frame_size; % number of data bits (overall rate is 1/2)
SNR = 10^(0.1*SNR_dB); % SNR per bit in linear scale
noise_var_1D = 2*2/(2*SNR); % 1D noise variance
%--------------------------------------------------------------------------
%    Generator polynomial of the inner encoder
gen_poly_inner = ldiv2([1 0 1],[1 1 1],2*num_bit); % using long division method

%    Generator polynomial of the outer encoder
gen_poly_outer = ldiv2([1 0 1],[1 1 1],num_bit); % using long division method

%  Interleaver and deinterleaver mapping of the SCCC 
intr_map = randperm(2*num_bit);
deintr_map = deintrlv((1:2*num_bit),intr_map);

%--------------------------------------------------------------------------
C_Ber = 0; % channel errors
tic()
%--------------------------------------------------------------------------
for frame_cnt = 1:sim_runs
%                           TRANSMITTER
%Source
a = randi([0 1],1,num_bit); % data

% SCCC encoder
% Outer encoder
b = zeros(1,2*num_bit); % outer encoder output initialization
b(1:2:end) = a; % systematic bit
temp1 = mod(conv(gen_poly_outer,a),2); % linear convolution with the generator polynomial
b(2:2:end) = temp1(1:num_bit); % parity bit

% interleaver
c = b(intr_map);

% Inner encoder
d = zeros(1,2*frame_size); % inner encoder output initialization
d(1:2:end) = c; % systematic bit
temp2 = mod(conv(gen_poly_inner,c),2); % linear convolution with the generator polynomial
d(2:2:end) = temp2(1:frame_size); % parity bit

% QPSK mapping (according to the set partitioning principles)
mod_sig = 1-2*d(1:2:end) + 1i*(1-2*d(2:2:end));

%--------------------------------------------------------------------------
%                            CHANNEL   
% AWGN
white_noise = sqrt(noise_var_1D)*randn(1,frame_size)+1i*sqrt(noise_var_1D)*randn(1,frame_size); 
Chan_Op = mod_sig + white_noise; % Chan_Op stands for channel output
%--------------------------------------------------------------------------
%                          RECEIVER 

% Branch metrices for the inner BCJR
QPSK_SYM = zeros(4,frame_size);
QPSK_SYM(1,:) = (1+1i)*ones(1,frame_size);
QPSK_SYM(2,:) = (1-1i)*ones(1,frame_size);
QPSK_SYM(3,:) = (-1+1i)*ones(1,frame_size);
QPSK_SYM(4,:) = (-1-1i)*ones(1,frame_size);

Dist = zeros(4,frame_size);
 Dist(1,:)=abs(Chan_Op-QPSK_SYM(1,:)).^2;
 Dist(2,:)=abs(Chan_Op-QPSK_SYM(2,:)).^2;
 Dist(3,:)=abs(Chan_Op-QPSK_SYM(3,:)).^2;
 Dist(4,:)=abs(Chan_Op-QPSK_SYM(4,:)).^2;
 log_gamma = -Dist/(2*noise_var_1D); % log gamma
 
% a priori LLR for inner decoder for 1st iteration
LLR = zeros(1,frame_size);

% iterative logMAP decoding
LLR = log_BCJR_inner(LLR,log_gamma,frame_size); % outputs extrinsic information
LLR = log_BCJR_outer(LLR(deintr_map),num_bit); %1

LLR = log_BCJR_inner(LLR(intr_map),log_gamma,frame_size); 
LLR = log_BCJR_outer(LLR(deintr_map),num_bit); %2

LLR = log_BCJR_inner(LLR(intr_map),log_gamma,frame_size); 
LLR = log_BCJR_outer(LLR(deintr_map),num_bit); %3

LLR = log_BCJR_inner(LLR(intr_map),log_gamma,frame_size); 
LLR = log_BCJR_outer(LLR(deintr_map),num_bit); %4

LLR = log_BCJR_inner(LLR(intr_map),log_gamma,frame_size); 
LLR = log_BCJR_outer(LLR(deintr_map),num_bit); %5

LLR = log_BCJR_inner(LLR(intr_map),log_gamma,frame_size); 
LLR = log_BCJR_outer(LLR(deintr_map),num_bit); %6

LLR = log_BCJR_inner(LLR(intr_map),log_gamma,frame_size); 
LLR = log_BCJR_outer(LLR(deintr_map),num_bit); %7

LLR = log_BCJR_inner(LLR(intr_map),log_gamma,frame_size);
LLR = log_BCJR_outer_END(LLR(deintr_map),num_bit); % 8: outputs aposteriori information

% hard decision 
dec_data = LLR<0;
% 
 % Calculating total bit errors
C_Ber = C_Ber + nnz(dec_data-a); 
end

BER = C_Ber/(sim_runs*num_bit)
toc()

