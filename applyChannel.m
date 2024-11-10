function [lossSig, chanDelay] = applyChannel(sig, prm, channel, varargin)

        narginchk(3,4);
        numTx = prm.numTx;
        numRx = prm.numRx;
        if nargin>3
            % preSig, for data transmission
            preSig = varargin{1}; 
            sigPad = [preSig; zeros(prm.numPadZeros,numTx); ...
                      sig; zeros(prm.numPadZeros,numTx)];
        else
            % No preSig, for sounding
            preSig = []; 
            sigPad = [sig; zeros(prm.numPadZeros,numTx)];
        end

        %%
        % fadeSig = channel(sigPad);
        % chanDelay = info(channel).ChannelFilterDelay;
        % 
        % % Remove the preamble, if present
        % if ~isempty(preSig)
        %     fadeSig(1:(length(preSig)+prm.numPadZeros),:) = [];
        % end
        % lossSig = fadeSig;

        % lossSig = fadeSig/sqrt(db2pow(prm.spLoss));
        
       
        %%
        
        % 
        Ns = 100;          % Number of scatterers
        % Place scatterers randomly in a circle from the center
        posCtr = (prm.txPosition+prm.rxPosition)/2;
        radCtr = prm.rays{1}(1).PropagationDistance*0.45;
        scatBound = [posCtr(1)-radCtr posCtr(1)+radCtr; ...
            posCtr(2)-radCtr posCtr(2)+radCtr;     0 0];

        % Channel
        channel = phased.ScatteringMIMOChannel(...
            'TransmitArray',prm.txarray,...
            'ReceiveArray',prm.rxarray,...
            'PropagationSpeed',prm.cLight,...
            'CarrierFrequency',prm.fc,...
            'SampleRate',prm.chanSRate, ...
            'SimulateDirectPath',false, ...
            'ChannelResponseOutputPort',true, ...
            'TransmitArrayPosition',prm.txPosition',...
            'ReceiveArrayPosition',prm.rxPosition',...
            'NumScatterers',Ns, ...
            'ScattererPositionBoundary',scatBound, ...
            'SeedSource','Property');

        [fadeSig, ~, tau] = channel(sigPad);
        chanDelay = floor(min(tau)*prm.chanSRate);

        % Remove the preamble, if present
        if ~isempty(preSig)
            fadeSig(1:(length(preSig)+prm.numPadZeros),:) = [];
        end

        % Path loss is included in channel
        lossSig = fadeSig;


end

