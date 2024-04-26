function [RX_startbyte, RX_endpointID, RX_content] = receiveMessage(obj)
    maxRetries = 10;  % Define the maximum number of retries for receiving the correct start byte
    retries = 0;      % Initialize the retry counter

    while retries < maxRetries
        % Receive start byte
        RX_startbyte = fread(obj.hSerialPort, 1, 'uint8');
        
        % Check if the start byte is as expected
        if (RX_startbyte == obj.cStartbyteMessage) || (RX_startbyte == obj.cStartbyteStatus)
            break;  % Correct start byte found, break the loop to proceed
        else
            % Unexpected start byte, increment retry counter and flush input buffer
            retries = retries + 1;
            flushinput(obj.hSerialPort);
            disp(['[RadarSystem.receiveMessage] Warning: Unexpected start byte received. Retrying... (' num2str(retries) '/' num2str(maxRetries) ')']);
        end
    end

    % Check if max retries have been exceeded
    if retries >= maxRetries
        error('[RadarSystem.receiveMessage] Error: Maximum retries exceeded while waiting for correct start byte.');
    end

    % If the correct start byte was found, proceed to read the endpoint ID
    RX_endpointID = fread(obj.hSerialPort, 1, 'uint8');

    % Proceed based on the start byte received
    if (RX_startbyte == obj.cStartbyteMessage)
        % Receive payload size
        RX_payloadSize = fread(obj.hSerialPort, 1, 'uint16');
        % Receive payload
        RX_content = cast(fread(obj.hSerialPort, RX_payloadSize, 'uint8'),'uint8')';
        % Receive message end sequence
        RX_payloadEODB = fread(obj.hSerialPort, 1, 'uint16');

        % Check if message was received correctly
        if (RX_payloadEODB ~= obj.cEndOfMessage)
            disp('[RadarSystem.receiveMessage] Error: Bad message end sequence received');
            RX_content = [];  % Clear content in case of error
        end
        
    elseif (RX_startbyte == obj.cStartbyteStatus)
        % Read the status content if the start byte indicates a status message
        RX_content = fread(obj.hSerialPort, 1, 'uint16');
        
    else
        disp('[RadarSystem.receiveMessage] Error: Bad message start byte received');
        RX_content = [];  % Clear content in case of error
    end
end
