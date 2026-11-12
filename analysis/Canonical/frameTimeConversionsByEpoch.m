function timeFrameRangeEpochs = frameTimeConversionsByEpoch(framesIn, timesIn, sz_t, frameRate, epochNames)
    %takes a range of frame numbers corresponding to experimental epochs.
    %ex: framesIn = [1 279; 280 510];  
    % (2 values each row, [startFrame endFrame]
    % Leave Empty [], if calculating from times

    %*** or takes a time range, corresponding to experimental epoch
    %ex: timeIn = [0 5; 5 10; 10 30]. (minutes)
    % Leave Empty [], if calculating from frames

    %sz_t is the number of timePoints in an image stack
    %frameRate is the timeInterval between stacks, given in seconds

    %takes optional text to describe the epochs
    %ex: epochNames = {'Baseline'; 'ExperimentalCondition'}

    %define defaults
    frameMin = 1;
    frameMax = sz_t;
    defaultFrameRange = [frameMin frameMax];
    %time step (in minutes)
    timeStep = frameRate / 60;
    %last time (in minutes)
    timeMax = timeStep*frameMax;
    defaultTimeRange = {[timeStep timeMax]};

    %check if TimeRange was input.  If so, use that to compute frameRange
    if ~isempty(timesIn) %time is given (priority)
        timeRange = timesIn;
        frameRange = frameRangeFromTimeRange(timesIn, timeStep, frameMin, frameMax);
    else %no time input
        if ~isempty(framesIn) %time was input
            timeRange = timeRangeFromFrameRange(framesIn, timeStep, timeStep, timeMax);
            frameRange = framesIn;
        else %no time input, no frame input, use defaults
            timeRange = defaultTimeRange;
            frameRange = defaultFrameRange;
        end
    end
    %check if epochNames are given
    nRanges = size(frameRange, 1);
    nEpochNames = length(epochNames);
    if ~(nRanges == nEpochNames)
        for i=1:nRanges
            epochNames{i} = strcat('Epoch', num2str(i)); 
        end
    end
    %outputCell
    timeFrameRangeEpochs.timeRange = timeRange;
    timeFrameRangeEpochs.frameRange = frameRange;
    timeFrameRangeEpochs.epochNames = epochNames;
   
end




function frameRange = frameRangeFromTimeRange(timeRange, timeStep, frameMin, frameMax)

    numEpochs = size(timeRange,1);
    frameRange = zeros(size(timeRange));
    for epoch = 1:numEpochs
        %this epoch's times
        thisT1 = timeRange(epoch, 1);
        thisT2 = timeRange(epoch, 2);
        %convert to framenumber
        if thisT1==0
            thisF1 = frameMin;
        else
            thisF1 = ceil(thisT1/timeStep);
        end
        if isinf(thisT2)
            thisF2 = frameMax;
        else
            thisF2 = floor(thisT2/timeStep);
        end
        %check limits
        thisF1 = max(thisF1, frameMin);
        thisF2 = min(thisF2, frameMax);
        %check overlap from previous
        if epoch>1
            previousF2 = frameRange(epoch, 2);
            if thisF1==previousF2
                thisF1 = thisF1 + 1;
            end
        end
        %store it
        frameRange(epoch, 1) = thisF1;
        frameRange(epoch, 2) = thisF2;
    end
end


function timeRange = timeRangeFromFrameRange(frameRange, timeStep, timeMin, timeMax)
    numEpochs = size(frameRange,1);
    timeRange = zeros(size(frameRange));
    for epoch = 1:numEpochs
        %this epoch's frames
        thisF1 = frameRange(epoch, 1);
        thisF2 = frameRange(epoch, 2);
        %convert to framenumber
        if thisF1==0
            thisT1 = timeMin;
        else
            thisT1 = thisF1*timeStep;
        end
        if isinf(thisF2)
            thisT2 = timeMax;
        else
            thisT2 = thisF2*timeStep;
        end
        %check limits
        thisT1 = max(thisT1, timeMin);
        thisT2 = min(thisT2, timeMax);
        %check overlap from previous
        if epoch>1
            previousT2 = timeRange(epoch, 2);
            if thisT1==previousT2
                thisT1 = thisT1 + timeStep;
            end
        end
        %store it
        timeRange(epoch, 1) = thisT1;
        timeRange(epoch, 2) = thisT2;
    end
end
