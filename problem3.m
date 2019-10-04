%**************************************************************************
%
% ASA Student Challenge Problem 2019
%
%   Script loads the time series hydrophone data and compares the measured
%   time frequency data with the calulated shifts.
%
%        Scott Schoen Jr | Georgia Tech | 30 September 2019
%
%**************************************************************************

clear all
close all
clc

% Problem parameters
d = 90; % [m]
h = 150; % [m]
v = 94; % [m/s]
f0 = 68.36; % [Hz]

c1 = 343; % [m/s]
c2 = 1520; % [m/s]

% Load in hydrophone data
dataFile = '../media/hydrophoneSignal.wav';
[hydSig, Fs] = audioread(dataFile);

% Downsample for efficiency
fMax = 30E3; % [Hz]
N = length( hydSig );
fRatio = Fs./fMax;
indsToKeep = round( 1 : fRatio : N );
hydSig = hydSig( indsToKeep );

% Set time of plane overhead
tOffset = 54.6; % [s]
measuredShift = 4.12; % [Hz]

% Normalize
hydSigNorm = hydSig./max(abs(hydSig));

% Get time signal
N = length( hydSig );
dt = 1./fMax;
tVec = 0 : dt : (N - 1).*dt;


% First sweep through and find the velocity to match this shift
vVec = linspace( 94.5, 96.5, 50 );
deltaf = 0.*vVec;

for vCount = 1 : length( vVec )
    
    % Define plane position at each point
    t = linspace( -4, 4, 100 );
    x = -vVec(vCount).*t;
    
    % Initialize
    thetac = asin( c1./c2 );
    theta0 = 0.99.*[ -thetac, thetac ];
    theta1 = 0.*t;
    
    % Find theta1 as function of time
    for tCount = 1 : length(t)
        
        % Current plane position
        D = x(tCount);
        
        % Function to be solved for theta1 at each point
        thetaVec = linspace( theta0(1), theta0(2), 100 );
        zeroFunction = @(thetaVar) ...
            D - (h.*tan(thetaVar) + d.*tan(asin((c2./c1).*sin(thetaVar))));
        theta1(tCount) = fzero( zeroFunction, theta0 );
        
    end
    
    % Compute theta2 depression angles
    theta2 = asin( (c2./c1).*sin(theta1) );
    gamma = pi./2 - theta1;
    xi = pi./2 - theta2;
    
    % Compute the instantaneous phase
    fw = f0.*( 1 - (v./c1).*cos( gamma ) ).^(-1); % At the surface
    
    % Get frequency range
    deltaf( vCount ) = max(fw) - min(fw);
    
end

% Get the optimal velocity
[~, vInd] = min( abs( deltaf./2 - measuredShift ) );

% Compute shifts for that velocity
v = vVec(vInd);
x = -vVec(vInd).*t;

% Initialize
thetac = asin( c1./c2 );
theta0 = 0.99.*[ -thetac, thetac ];
theta1 = 0.*t;

% Find theta1 as function of time
for tCount = 1 : length(t)
    
    % Current plane position
    D = x(tCount);
    
    % Function to be solved for theta1 at each point
    thetaVec = linspace( theta0(1), theta0(2), 100 );
    zeroFunction = @(thetaVar) ...
        D - (h.*tan(thetaVar) + d.*tan( asin( (c2./c1).*sin(thetaVar) ) ) );
    theta1(tCount) = fzero( zeroFunction, theta0 );
    
end

% Compute theta2 depression angles
theta2 = asin( (c2./c1).*sin(theta1) );
gamma = pi./2 - theta1;
xi = pi./2 - theta2;

% Compute the instantaneous phase
fw = f0.*( 1 - (v./c1).*cos( gamma ) ).^(-1); % At the surface


%% Plot spectrogram and overlay the predicted frequency
figure();
hold all;

[s, fVec, tVec] = spectrogram(hydSigNorm, 2.^(15), 2.^(14), [], fMax, 'yaxis');
[T, F] = meshgrid( tVec, fVec );

pcolor( T - tOffset, F, 20.*log10(abs(s)./max(abs(s(:)))) );
colormap( flipud( colormap( gray ) ) );
shading interp;

plot( t, fw, 'k', 'LineWidth', 5 );

ylim( [55, 85] );
xlim(3.75.*[ -1, 1]);
caxis([-24, 0]);

xlabel( 'Time [s]', 'FontSize', 24 );
ylabel( 'Frequency [Hz]', 'FontSize', 24 );

set( gca, 'XTick', -5:5 );

cbh = colorbar;
cbh.Location = 'North';
cbh.Position = [0.675, 0.8, 0.17, 0.02];
cbh.Ticks = [];


grid on;
box off;

