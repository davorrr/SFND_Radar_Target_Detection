clear all
clc;

%% Radar Specifications 
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Frequency of operation = 77GHz
% Max Range = 200m
range_max = 200;
% Range Resolution = 1 m
ds = 1;
% Max Velocity = 100 m/s
%%%%%%%%%%%%%%%%%%%%%%%%%%%

%speed of light = 3e8
c = 3e8;
%% User Defined Range and Velocity of target
% *%TODO* :
% define the target's initial position and velocity. Note : Velocity
% remains contant

target_position = 80;
velocity = 20;
 


%% FMCW Waveform Generation

% *%TODO* :
% Design the FMCW waveform by giving the specs of each of its parameters.
% Calculate the Bandwidth (B), Chirp Time (Tchirp) and Slope (slope) of the FMCW
% chirp using the requirements above.
B = c / 2* ds;
Tchirp = 5.5 * (range_max*2/c);
slope = B / Tchirp;
disp(slope);

% Operating carrier frequency of Radar 
fc= 77e9;             %carrier freq
                                                          
% The number of chirps in one sequence. Its ideal to have 2^ value for the ease of running the FFT
%for Doppler Estimation. 
Nd=128;                   % #of doppler cells OR #of sent periods % number of chirps

% The number of samples on each chirp. 
Nr=1024;                  %for length of time OR # of range cells

% Timestamp for running the displacement scenario for every sample on each
% chirp
t=linspace(0,Nd*Tchirp,Nr*Nd); %total time for samples

% Creating the vectors for Tx, Rx and Mix based on the total samples input.
Tx=zeros(1,length(t)); %transmitted signal
Rx=zeros(1,length(t)); %received signal
Mix = zeros(1,length(t)); %beat signal

% Similar vectors for range_covered and time delay.
r_t=zeros(1,length(t));
td=zeros(1,length(t));

%% Signal generation and Moving Target simulation
% Running the radar scenario over the time. 

for i=1:length(t)         
        
    % *%TODO* :
    % For each time stamp update the Range of the Target for constant velocity. 
    r_t(i) = target_position + t(i)*velocity;
    td(i) = 2 * r_t(i)/c;

    % *%TODO* :
    % For each time sample we need update the transmitted and
    % received signal. 
    Tx(i) = cos(2 * pi * (fc * t(i) + 0.5 * slope * t(i)^2));
    Rx(i) = cos(2 * pi * (fc * (t(i) - td(i)) + 0.5 * slope * (t(i) - td(i))^2));
  
    % *%TODO* :
    % Now by mixing the Transmit and Receive generate the beat signal
    % This is done by element wise matrix multiplication of Transmit and
    % Receiver Signal
    Mix(i) = Tx(i) .* Rx(i);
    
end

%% RANGE MEASUREMENT

% *%TODO* :
% reshape the vector into Nr*Nd array. Nr and Nd here would also define the size of
% Range and Doppler FFT respectively.
Mix = reshape(Mix, [Nr, Nd]);

% *%TODO* :
%run the FFT on the beat signal along the range bins dimension (Nr) and
%normalize.
% TODO : Compute the Fourier transform of the signal. 
signal_fft = fft(Mix);

% *%TODO* :
% Take the absolute value of FFT output
L = Tchirp * B;
signal_fft = abs(signal_fft/L);

% *%TODO* :
% Output of FFT is double sided signal, but we are interested in only one side of the spectrum.
signal_fft = signal_fft(1:L/2+1);

%plotting the range
f = B * (0:(L/2))/L;
R = (c * Tchirp * f) / (2 * B);
figure ('Name','Range from First FFT')
plot(R, signal_fft);
xlabel('Doppler velocity [m/s]')
xlabel('Range [m]')
ylabel('Doppler velocity [m/s]')
axis ([0 200 0 0.3]);

%TODO* :
% plot FFT output 
figure ('Name','FFT')
plot(f, signal_fft); 




%% RANGE DOPPLER RESPONSE

% Range Doppler Map Generation.

% The output of the 2D FFT is an image that has reponse in the range and
% doppler FFT bins. So, it is important to convert the axis from bin sizes
% to range and doppler based on their Max values.
Mix=reshape(Mix,[Nr,Nd]);

% 2D FFT using the FFT size for both dimensions.
sig_fft2 = fft2(Mix,Nr,Nd);

% Taking just one side of signal from Range dimension.
sig_fft2 = sig_fft2(1:Nr/2,1:Nd);
sig_fft2 = fftshift (sig_fft2);
RDM = abs(sig_fft2);
RDM = 10*log10(RDM) ;

% Use the surf function to plot the output of 2DFFT and to show axis in both
% dimensions
doppler_axis = linspace(-100,100,Nd);
range_axis = linspace(-200,200,Nr/2)*((Nr/2)/400);
figure,surf(doppler_axis,range_axis,RDM)
title('Range Doppler Response')
xlabel('Doppler velocity [m/s]')
ylabel('Range [m]')
zlabel('RxPwr [dBW]')

%% CFAR implementation

% Slide Window through the complete Range Doppler Map

% *%TODO* :
% Create a vector to store noise_level for each iteration on training cells
signal_cfar_2D = zeros(Nr/2, Nd);
noise_level = zeros(1,1);

% *%TODO* :
% design a loop such that it slides the CUT across range doppler map by
% giving margins at the edges for Training and Guard Cells.
% For every iteration sum the signal level within all the training
% cells. To sum convert the value from logarithmic to linear using db2pow
% function. Average the summed values for all of the training
% cells used. After averaging convert it back to logarithimic using pow2db.
% Further add the offset to it to determine the threshold. Next, compare the
% signal under CUT with this threshold. If the CUT level > threshold assign
% it a value of 1, else equate it to 0.

% Use RDM[x,y] as the matrix from the output of 2D FFT for implementing CFAR

% The number of Training Cells in both the dimensions.
Tr = 14; % Training cells in Range dimension
Td = 8; % Training cells in Doppler dimension

% The number of Guard Cells in both dimensions around the Cell under 
% test (CUT) for accurate estimation
Gr = 8; % Guard cells in Range dimension
Gd = 4; % Guard cells in Doppler dimension

% The offset the threshold by SNR value in dB
offset = 7;

% Total number of training cells
num_cells = (2*Tr + 2*Gr + 1)*(2*Td + 2*Gd + 1) - (2*Gr + 1)*(2*Gd + 1);

% Nr - number of Range cells 
% Nd - number of Doppler cells
for i = 1:(Nr/2 - (2*Gr + 2*Tr))  % looping across cells  in range dimension (X)
    for j = 1:(Nd - (2*Gd + 2*Td)) % for each cell in range dimension looping across cells in doppler dimension (Y)

        % For every iteration sum the signal level within all the training
        % cells. To sum convert the value from logarithmic to linear using 
        % db2pow function.
        signal_1 = sum(db2pow(RDM(i : i + 2*Tr + 2*Gr, j : j + 2*Td + 2*Gd)), "all");
        signal_2 = sum(db2pow(RDM(i + Tr : i + Tr + 2*Gr, j + Td : j + Td + 2*Gd)), "all");
        signal_level = signal_1 - signal_2;

        % Average the summed values for all of the training cells used. 
        % After averaging convert it back to logarithmic using pow2db. 
        % Measure and average the noise across all the training cells. This
        % gives the threshold. Add the offset (if in signal strength in dB)
        % to the threshold to keep the false alarm to the minimum
        threshold = signal_level / num_cells;
        threshold = pow2db(threshold) + offset;
        threshold = db2pow(threshold);

        % Determine the signal level at the Cell Under Test.
        signal = db2pow(RDM(i + Tr + Gr, j + Td + Gd));

        % If the CUT signal level is greater than the Threshold, assign a 
        % value of 1, else equate it to zero.
        if(signal <= threshold)
            signal = 0;
        else
            signal = 1;
        end
     
        signal_cfar_2D(i + Tr + Gr, j + Td + Gd) = signal;
    end
end

% *%TODO* :
% display the CFAR output using the Surf function like we did for Range
% Doppler Response output.
figure,surf(doppler_axis,range_axis,signal_cfar_2D);
colorbar;
title('2D CA-CFAR')
xlabel('Doppler velocity [m/s]')
ylabel('Range [m]')
zlabel('Detection')

