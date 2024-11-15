function prm = array_setup(prm,varargin)
% Return whether Tx and Rx employ URA or not, with other information
%
        narginchk(1,2);
        if nargin==2
            chkFlag = varargin{1};
        else
            chkFlag = false;
        end
        
        if chkFlag
            % Parameter checks
            %   numSTS
            if (~ismember(prm.numSTS, [1 2 4 8 16 32 64]))
                error('The number of data streams must be one of [1 2 4 8 16 32 64].');
            end
            
            %   numSTS and numTx
            if (prm.numSTS > prm.numTx)
                error('The number of data streams must be less than or equal to the number of transmit antennas.');
            end
        
            if (floor(prm.numTx/prm.numSTS) ~= prm.numTx/prm.numSTS)
                error('The number of transmit antennas must be an integer multiple of the total number of data streams.');
            end
            
            %   numSTS and numRx
            if any(prm.numRx < prm.numSTSVec)
                error('The number of receive antennas per user must be greater than or equal to the number of data streams per user.');
            end
        end
        
        % Transmit end (Tx)
        expFactorTx = prm.numTx/prm.numSTS;   
        
        if expFactorTx==1 || prm.numSTS==1 % ULA
            % Use ULA if expFactorTx==1 or if prm.numSTS<2
            isTxURA = false;
            if prm.numTx==1
                % Avoid the phased.ULA error for numElements==1
                error('The number of transmit antennas must be greater than 1.');
            end
        else % URA
            % Use URA if expFactorTx>1 and numSTS>1
            isTxURA = true;
        end
        
        % Receive end (Rx), per user
        expFactorRx = prm.numRx./prm.numSTSVec;
        
        isRxURA = zeros(prm.numUsers, 1, 'logical');
        for uIdx = 1:prm.numUsers
            if floor(expFactorRx(uIdx))~=expFactorRx(uIdx) || expFactorRx(uIdx)==1 ...
                    || prm.numSTSVec(uIdx)==1   % ULA
                % Use ULA if expFactorRx==1 or expFactorRx is non-integer
                isRxURA(uIdx) = false;
            else % URA
                % Use URA if expFactorRx>1 and an integer, and numSTS>1
                isRxURA(uIdx) = true;
            end
        end



        
        prm.isTxURA = isTxURA;         % Assign isTxURA to its value
        prm.expFactorTx = expFactorTx; % Assign expFactorTx to its value
        prm.isRxURA = isRxURA;         % Assign isRxURA to its value
        prm.expFactorRx = expFactorRx; % Assign expFactorRx to its value
        prm.txPosition = [23.726582,90.389145, 23]; % lat, lon, height
        prm.txArrayOrientation = [-90 0].';       % azimuth (0 deg is East, 90 deg is North) and elevation (positive points upwards) in deg
        prm.rxPosition = [23.725462,90.387456, 15];   % lat, lon
        prm.rxArrayOrientation = [90 0].';      % azimuth (0 deg is East, 90 deg is North) and elevation (positive points upwards)  in deg
        prm.reflectionsOrder = 1;                 % number of reflections for ray tracing analysis (0 for LOS)
        %% Phased array antennar setup
        if prm.isTxURA % URA
            prm.txarray = phased.URA([prm.expFactorTx,prm.numSTS], ...
                [0.5 0.5]*prm.lambda,'Element', ...
                phased.IsotropicAntennaElement('BackBaffled',true));
        else % ULA
            prm.txarray = phased.ULA('Element', ...
                phased.IsotropicAntennaElement('BackBaffled', true),...
                'NumElements',prm.numTx,'ElementSpacing',0.5*prm.lambda);
        end         
        if prm.isRxURA % URA
            prm.rxarray = phased.URA([expFactorRx,prm.numSTS], ...
                [0.5 0.5]*prm.lambda,'Element', ...
                phased.IsotropicAntennaElement('BackBaffled',true));
        else % ULA
            prm.rxarray = phased.ULA('Element',phased.IsotropicAntennaElement, ...
                'NumElements',prm.numRx,'ElementSpacing',0.5*prm.lambda);
        end


        %% Locate the BS and UE on the site
        prm.txSite = txsite("Name","Base station", ...
            "Latitude",prm.txPosition(1),"Longitude",prm.txPosition(2),...
            "AntennaAngle",prm.txArrayOrientation(1:2),...
            "AntennaHeight",prm.txPosition(3),...  % in m
            "Antenna", prm.txarray,...
            "TransmitterFrequency",prm.fc);
        
        prm.rxSite = rxsite("Name","UE", ...
            "Latitude",prm.rxPosition(1),"Longitude",prm.rxPosition(2),...
            "AntennaHeight",prm.rxPosition(3),... % in m
            "Antenna", prm.rxarray,...
            "AntennaAngle",prm.rxArrayOrientation(1:2));
        %% Ray tracing analysis
        pm = propagationModel("raytracing","Method","sbr","MaxNumReflections",prm.reflectionsOrder);
        prm.rays = raytrace(prm.txSite,prm.rxSite,pm,"Type","pathloss");   

end
