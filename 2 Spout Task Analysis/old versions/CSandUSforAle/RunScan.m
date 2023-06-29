function RunScan()
S.arduinoPort = 'COM12';
S.wheelCirc = 46; % Circumference of the wheel (cm)
S.wheelTicks = 1024; % Number of ticks recorded by the wheel per rotation.
S.wheelMult = S.wheelCirc/S.wheelTicks; % Distance per tick (cm)
%Attempt to connect to the Arduino run detector.
ardSerial = serial(S.arduinoPort);
ardSerial.InputBufferSize = 16000; % read up to 2kb
set(ardSerial,'DataBits',8);
set(ardSerial,'StopBits',1);
set(ardSerial,'BaudRate',115200);
set(ardSerial,'Parity','none');
fopen(ardSerial);
autoCleanup = onCleanup(@()fclose(ardSerial));
%% Create a figure window to monitor the live data
Tmax = 3; % Total time for data collection (s)
%figure,
%grid on,
%xlabel ('Time (s)'), ylabel('Data (8-bit)'),
%axis([0 Tmax+1 -10 300]),
%% Read and plot the data from Arduino
Ts = 1; % Sampling time (s)
i = 0;
data = [];
t = 0;
tic % Start timer
while toc <= Tmax
	i = i + 1;
	%% Read buffer data
	%data(i) = fread(ardSerial,1000);
	dataSize = get(ardSerial,'BytesAvailable');
	if (dataSize>0)
		inASCII = fread(ardSerial,dataSize);
		inString = transpose(char(inASCII));
		inString = strsplit(inString);
		inString = inString(1:end-1); % Remove the last newline from the array
		for k=1:length(inString)
			disp(inString(k));
			inString{k} = str2double(strsplit(inString{k},'A'));
		end
		assignin('base', 'inString', inString);
		disp(inData);
		disp('\n')
		disp(i)
	end
	%% Read time stamp
	% If reading faster than sampling rate, force sampling time.
	% If reading slower than sampling rate, nothing can be done. Consider
	% decreasing the set sampling time Ts
	t(i) = toc;
	if i > 1
		T = toc - t(i-1);
		while T < Ts
			T = toc - t(i-1);
		end
	end
	t(i) = toc;
	%% Plot live data
	if i > 1
		line([t(i-1) t(i)],[data(i-1) data(i)])
		drawnow
	end
end
end