

close all;
clc;
s = rng(61);        % Set RNG state for repeatability

%% System Parameters
%
% Define parameters for the system. These parameters can be modified to
% explore their impact on the system.

% Single-user system with multiple streams
prm.numUsers = 1;            % Number of users
prm.numSTS = 16;             % Number of independent data streams, 4/8/16/32/64
prm.numTx = 32;              % Number of transmit antennas 
prm.numRx = 16;              % Number of receive antennas 
prm.bitsPerSubCarrier = 6;   % 2: QPSK, 4: 16QAM, 6: 64QAM, 8: 256QAM
prm.numDataSymbols = 10;     % Number of OFDM data symbols

prm.fc = 4e9;                   % 4 GHz system
prm.chanSRate = 100e6;          % Channel sampling rate, 100 Msps
prm.NFig = 5;                   % Noise figure, dB
prm.enSteering = false;          % Enable/disable steering
prm.mobileAngle = [-90; 0];

% Parameters to define the OFDM modulation employed for the system are
% specified below.

prm.FFTLength = 256; 
prm.CyclicPrefixLength = 64; 
prm.numCarriers = 234; 
prm.NumGuardBandCarriers = [7 6];
prm.PilotCarrierIndices = [26 54 90 118 140 168 204 232];
nonDataIdx = [(1:prm.NumGuardBandCarriers(1))'; prm.FFTLength/2+1; ...
              (prm.FFTLength-prm.NumGuardBandCarriers(2)+1:prm.FFTLength)'; ...
              prm.PilotCarrierIndices.';];
prm.CarriersLocations = setdiff((1:prm.FFTLength)',sort(nonDataIdx));

numTx = prm.numTx;
numRx = prm.numRx;
numSTS = prm.numSTS;
prm.numFrmBits = numSTS*prm.numDataSymbols*prm.numCarriers* ...
                 prm.bitsPerSubCarrier*1/3-6; % Account for termination bits

prm.modMode = 2^prm.bitsPerSubCarrier; % Modulation order

% Account for channel filter delay
prm.numPadZeros = 3*(prm.FFTLength+prm.CyclicPrefixLength); 

% Get transmit and receive array information
prm.numSTSVec = numSTS;
prm.cLight = physconst('LightSpeed');
prm.lambda = prm.cLight/prm.fc;
gainFactor = 1;

%% Array Placement
prm = array_setup(prm);

%% Import and Visualize 3-D Environment with Buildings for Ray Tracing
if exist('viewer','var') && isvalid(viewer) % viewer handle exists and viewer window is open
    clearMap(viewer);
else
    viewer = siteviewer("Basemap","openstreetmap","Buildings","dhaka.osm");    
end

show(prm.txSite);
show(prm.rxSite);
plot(prm.rays{1});
%%

prm.avgPathGains  = -[prm.rays{1}.PathLoss];                                    % Average path gains of each ray
% prm.spLoss = -mean(prm.avgPathGains);
prm.spLoss = -prm.avgPathGains(1);
prm.pathAoDs = [prm.rays{1}.AngleOfDeparture];                                  % AoD of each ray
prm.pathAoAs = [prm.rays{1}.AngleOfArrival];                                    % AoA of each ray
prm.isLOS = any([prm.rays{1}.LineOfSight]);  
%% Channel Sounding


% Generate the preamble signal
preambleSigSTS = helperGenPreamble(prm);
%   repeat over numTx
preambleSig = zeros(size(preambleSigSTS,1),numTx);
for i = 1:numSTS
    preambleSig(:,(i-1)*prm.expFactorTx+(1:prm.expFactorTx)) = ...
        repmat(preambleSigSTS(:,i),1,prm.expFactorTx);
end
% get the channel 
[prm, channel] = get_channel(prm);
% Transmit preamble over channel
[rxPreSig,chanDelay] = applyChannel(preambleSig, prm, channel);


% Front-end amplifier gain and thermal noise
rxPreAmp = phased.ReceiverPreamp( ...
    'Gain',gainFactor*prm.spLoss, ... % account for path loss
    'NoiseFigure',prm.NFig, ...
    'ReferenceTemperature',290, ...
    'SampleRate',prm.chanSRate);
rxPreSigAmp = rxPreAmp(rxPreSig);
rxPreSigAmp = rxPreSigAmp * ...         % scale power
    (sqrt(prm.FFTLength-sum(prm.NumGuardBandCarriers)-1)/(prm.FFTLength));  

% OFDM Demodulation
demodulatorOFDM = comm.OFDMDemodulator( ...
     'FFTLength',prm.FFTLength, ...
     'NumGuardBandCarriers',prm.NumGuardBandCarriers.', ...
     'RemoveDCCarrier',true, ...
     'PilotOutputPort',true, ...
     'PilotCarrierIndices',prm.PilotCarrierIndices.', ...
     'CyclicPrefixLength',prm.CyclicPrefixLength, ...
     'NumSymbols',numSTS, ... % preamble symbols alone
     'NumReceiveAntennas',numRx);

rxOFDM = demodulatorOFDM( ...
    rxPreSigAmp(chanDelay+1:end-(prm.numPadZeros-chanDelay),:));

% Channel estimation from preamble
%       numCarr, numSTS, numRx
hD = helperMIMOChannelEstimate(rxOFDM(:,1:numSTS,:),prm); 

% Calculate the feedback weights
v = diagbfweights(hD);


%% Data Transmission 
%


% Convolutional encoder
encoder = comm.ConvolutionalEncoder( ...
    'TrellisStructure',poly2trellis(7,[133 171 165]), ...
    'TerminationMethod','Terminated');

% Generate mapped symbols from bits
txBits = randi([0, 1],prm.numFrmBits,1);
% txBits = ones(prm.numFrmBits,1);
encodedBits = encoder(txBits);

% Bits to QAM symbol mapping
mappedSym = qammod(encodedBits,prm.modMode,'InputType','Bit', ...
    'UnitAveragePower',true);

% Map to layers: per symbol, per data stream
gridData = reshape(mappedSym,prm.numCarriers,prm.numDataSymbols,numSTS);

% Apply normalized precoding weights to the subcarriers, assuming perfect feedback
P = complex(zeros(numSTS,numSTS,prm.numCarriers));
for carrIdx = 1:prm.numCarriers
    Q = squeeze(v(carrIdx,:,:));
    normQ = Q * sqrt(numTx)/norm(Q,'fro');  % normalize    
    P(:,:,carrIdx) = normQ;
end
preData = ofdmPrecode(gridData,P);

% OFDM modulation of the data
modulatorOFDM = comm.OFDMModulator( ...
    'FFTLength',prm.FFTLength,...
    'NumGuardBandCarriers',prm.NumGuardBandCarriers.',...
    'InsertDCNull',true, ...
    'PilotInputPort',true,...
    'PilotCarrierIndices',prm.PilotCarrierIndices.',...
    'CyclicPrefixLength',prm.CyclicPrefixLength,...
    'NumSymbols',prm.numDataSymbols,...
    'NumTransmitAntennas',numSTS);

% Multi-antenna pilots
pilots = helperGenPilots(prm.numDataSymbols,numSTS);

txOFDM = modulatorOFDM(preData,pilots);
txOFDM = txOFDM * (prm.FFTLength/ ...
    sqrt(prm.FFTLength-sum(prm.NumGuardBandCarriers)-1)); % scale power

% Generate preamble with the feedback weights and prepend to data
preambleSigD = helperGenPreamble(prm,v);
txSigSTS = [preambleSigD;txOFDM];

% Repeat over numTx
txSig = zeros(size(txSigSTS,1),numTx);
for i = 1:numSTS
    txSig(:,(i-1)*prm.expFactorTx+(1:prm.expFactorTx)) = ...
        repmat(txSigSTS(:,i),1,prm.expFactorTx);
end

%%


%% Transmit Beam Steering 

% Gain per antenna element 
amplifier = phased.Transmitter('PeakPower',1/numTx,'Gain',0);

% Amplify to achieve peak transmit power for each element
for n = 1:numTx
    txSig(:,n) = amplifier(txSig(:,n));
end


% For evaluating weights for steering  
SteerVecTx = phased.SteeringVector('SensorArray',prm.txarray, ...
    'PropagationSpeed',prm.cLight);

txSteerAngles = prm.txArrayOrientation-prm.pathAoDs;

% Generate weights for steered direction
wT = SteerVecTx(prm.fc,txSteerAngles(:,1));

% Radiate along the steered direction, without signal combining
radiatorTx = phased.Radiator('Sensor',prm.txarray, ...
    'WeightsInputPort',true, ...
    'PropagationSpeed',prm.cLight, ...
    'OperatingFrequency',prm.fc, ...
    'CombineRadiatedSignals',false);

if prm.enSteering
    txSteerSig = radiatorTx(txSig,repmat(prm.mobileAngle,1,numTx), ...
        conj(wT));
else
    txSteerSig = txSig;
end

% Visualize the array
h = figure('Position',figposition([10 55 22 35]),'MenuBar','none');
h.Name = 'Transmit Array Geometry';
viewArray(prm.txarray);

% Visualize the transmit pattern and steering
h = figure('Position',figposition([32 55 22 30]),'MenuBar','none');
h.Name = 'Transmit Array Response Pattern';
pattern(prm.txarray,prm.fc,'PropagationSpeed',prm.cLight,'Weights',wT);
h = figure('Position',figposition([54 55 22 35]),'MenuBar','none');
h.Name = 'Transmit Array Azimuth Pattern';
patternAzimuth(prm.txarray,prm.fc,'PropagationSpeed',prm.cLight,'Weights',wT);
if prm.isTxURA
    h = figure('Position',figposition([76 55 22 35]),'MenuBar','none');
    h.Name = 'Transmit Array Elevation Pattern';
    patternElevation(prm.txarray,prm.fc,'PropagationSpeed',prm.cLight, ...
        'Weights',wT);
end



%% Signal Propagation 
%

% Apply a spatially defined channel to the steered signal
[rxSig,chanDelay] = applyChannel(txSteerSig,prm,channel,preambleSig);

%% 
% The same channel is used for both sounding and data transmission, with
% the data transmission having a longer duration controlled by the number
% of data symbols parameter, |prm.numDataSymbols|.

%% Receive Beam Steering


rxPreAmp = phased.ReceiverPreamp( ...
    'Gain',gainFactor*prm.spLoss, ... % accounts for path loss
    'NoiseFigure',prm.NFig, ...
    'ReferenceTemperature',290, ...
    'SampleRate',prm.chanSRate);

% Front-end amplifier gain and thermal noise
rxSigAmp = rxPreAmp(rxSig);
rxSigAmp = rxSigAmp * ...           % scale power
    (sqrt(prm.FFTLength - sum(prm.NumGuardBandCarriers)-1)/(prm.FFTLength)); 

% For evaluating receive-side steering weights 
SteerVecRx = phased.SteeringVector('SensorArray',prm.rxarray, ...
    'PropagationSpeed',prm.cLight);

% Generate weights for steered direction towards mobile
% rxSteerAng = prm.pathAoAs(:, 1);
rxSteerAngles = prm.rxArrayOrientation-prm.pathAoAs;
wR = SteerVecRx(prm.fc,rxSteerAngles(:, 1));

% Steer along the mobile receive direction
if prm.enSteering
    rxSteerSig = rxSigAmp.*(wR');
else
    rxSteerSig = rxSigAmp;
end

% Visualize the array
h = figure('Position',figposition([10 20 22 35]),'MenuBar','none');
h.Name = 'Receive Array Geometry';
viewArray(prm.rxarray);

% Visualize the receive pattern and steering
h = figure('Position',figposition([32 20 22 30]));
h.Name = 'Receive Array Response Pattern';
pattern(prm.rxarray,prm.fc,'PropagationSpeed',prm.cLight,'Weights',wR);
h = figure('Position',figposition([54 20 22 35]),'MenuBar','none');
h.Name = 'Receive Array Azimuth Pattern';
patternAzimuth(prm.rxarray,prm.fc,'PropagationSpeed',prm.cLight,'Weights',wR);
if prm.isRxURA 
    figure('Position',figposition([76 20 22 35]),'MenuBar','none');
    h.Name = 'Receive Array Elevation Pattern';
    patternElevation(prm.rxarray,prm.fc,'PropagationSpeed',prm.cLight, ...
        'Weights',wR);
end

%%
% The receive antenna pattern mirrors the transmission steering.

%% Signal Recovery
%
% The receive antenna array passes the propagated signal to the receiver to
% recover the original information embedded in the signal. Similar to the
% transmitter, the receiver used in a MIMO-OFDM system contains many
% components, including OFDM demodulator, MIMO equalizer, QAM demodulator,
% and channel decoder.

demodulatorOFDM = comm.OFDMDemodulator( ...
     'FFTLength',prm.FFTLength, ...
     'NumGuardBandCarriers',prm.NumGuardBandCarriers.', ...
     'RemoveDCCarrier',true, ...
     'PilotOutputPort',true, ...
     'PilotCarrierIndices',prm.PilotCarrierIndices.', ...
     'CyclicPrefixLength',prm.CyclicPrefixLength, ...
     'NumSymbols',numSTS+prm.numDataSymbols, ... % preamble & data
     'NumReceiveAntennas',numRx);
  
% OFDM Demodulation
rxOFDM = demodulatorOFDM( ...
    rxSteerSig(chanDelay+1:end-(prm.numPadZeros-chanDelay),:));

% Channel estimation from the mapped preamble
hD = helperMIMOChannelEstimate(rxOFDM(:,1:numSTS,:),prm);

% MIMO Equalization
[rxEq,CSI] = ofdmEqualize(rxOFDM(:,numSTS+1:end,:),hD);

% Soft demodulation
scFact = ((prm.FFTLength-sum(prm.NumGuardBandCarriers)-1) ...
         /prm.FFTLength^2)/numTx;
nVar = noisepow(prm.chanSRate,prm.NFig,290)/scFact;
rxSymbs = rxEq(:)/sqrt(numTx);
rxLLRBits = qamdemod(rxSymbs,prm.modMode,'UnitAveragePower',true, ...
    'OutputType','approxllr','NoiseVariance',nVar);

% Apply CSI prior to decoding
rxLLRtmp = reshape(rxLLRBits,prm.bitsPerSubCarrier,[], ...
                   prm.numDataSymbols,numSTS);
csitmp = reshape(CSI,1,[],1,numSTS);
rxScaledLLR = rxLLRtmp.*csitmp;

% Soft-input channel decoding
decoder = comm.ViterbiDecoder(...
     'InputFormat','Unquantized', ...
     'TrellisStructure',poly2trellis(7, [133 171 165]), ...
     'TerminationMethod','Terminated', ...
     'OutputDataType','double');
rxDecoded = decoder(rxScaledLLR(:));

% Decoded received bits
rxBits = rxDecoded(1:prm.numFrmBits);

%% Received Signal Constellation 


% Display received constellation
constDiag = comm.ConstellationDiagram( ...
    'SamplesPerSymbol',1, ...
    'ShowReferenceConstellation',true, ...
    'ReferenceConstellation', ...
    qammod((0:prm.modMode-1)',prm.modMode,'UnitAveragePower',true), ...
    'ColorFading',false, ...
    'Position',figposition([20 20 35 40]), ...
    'Title','Equalized Symbols', ...
    'EnableMeasurements',true, ...
    'MeasurementInterval',length(rxSymbs));
constDiag(rxSymbs);

% Compute and display bit error rate
ber = comm.ErrorRate;
measures = ber(txBits,rxBits);
fprintf('BER = %.5f; No. of Bits = %d; No. of errors = %d\n', ...
    measures(1),measures(3),measures(2));

rng(s); % Restore RNG state
%%
% Plot UE radiation pattern
prm.rxSite.Antenna = clone(channel.ReceiveArray); % need a clone, otherwise setting the Taper weights would affect the channel array
prm.rxSite.Antenna.Taper = wR;
pattern(prm.rxSite,prm.fc,"Size",4);

% Plot BS radiation pattern
prm.txSite.Antenna = clone(channel.TransmitArray); % need a clone, otherwise setting the Taper weights would affect the channel array
prm.txSite.Antenna.Taper = wT;
pattern(prm.txSite,prm.fc,"Size",5);


%% Selected Bibliography
% # Perahia, Eldad, and Robert Stacey. Next Generation Wireless LANS:
% 802.11n and 802.11ac. Cambridge University Press, 2013.
% # IEEE(R) Std 802.11(TM)-2012 IEEE Standard for Information technology -
% Telecommunications and information exchange between systems - Local and
% metropolitan area networks - Specific requirements - Part 11: Wireless
% LAN Medium Access Control (MAC) and Physical Layer (PHY) Specifications.
% # 3GPP TS 36.213. "Physical layer procedures." 3rd Generation Partnership
% Project; Technical Specification Group Radio Access Network; Evolved
% Universal Terrestrial Radio Access (E-UTRA). URL: https://www.3gpp.org.
% # 3GPP TS 36.101. "User Equipment (UE) Radio Transmission and Reception."
% 3rd Generation Partnership Project; Technical Specification Group Radio
% Access Network; Evolved Universal Terrestrial Radio Access (E-UTRA). URL:
% https://www.3gpp.org.
% # Kyosti, Pekka, Juha Meinila, et al. WINNER II Channel Models.
% D1.1.2, V1.2. IST-4-027756 WINNER II, September 2007.
% # George Tsoulos, Ed., "MIMO System Technology for Wireless
% Communications", CRC Press, Boca Raton, FL, 2006.
