%**************************************************************************
%
% ASA Student Challenge Problem 2019
%
%   Problem 2 - Plots the  Doppler shifted frequency as a function of 
%   the bearing angle from the array for comparison with the provided
%   measurements.
%
%        Scott Schoen Jr | Georgia Tech | 30 September 2019
%
%**************************************************************************

clear all
close all
clc

dataFile = '../media/freqVsBearing.txt';

% Load in data
allData = importdata( dataFile );
fVec = allData( :, 1 );
bearingVec = allData(:, 2 );

% Problem parameters
d = 20; % [m]
h = 151; % [m]
v = 125; % [m/s]
f0 = 68.3; % [Hz]

c1 = 343; % [m/s]
c2 = 1520; % [m/s]

% Define plane position at each point
t = linspace( -1.4, 1.4, 100 );
x = -v.*t;

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


% Plot bearing vs frequency
figure();
hold all;

plot( 180.*xi./pi, fw, 'k', 'LineWidth', 3 ); 
plot( bearingVec, fVec, 'ko', 'MarkerFaceColor', 'k' );

set( gca, 'XTick', 0:30:180 );

ylabel( 'Frequency [Hz]', 'FontSize', 24 );
xlabel( 'Bearing Angle $\xi$ [deg]', 'FontSize', 24 );

ylim([60, 80]);
grid on;
box off;

lh = legend( '~~Prediction [Eq. (2.1)]', '~~Measurement', ...
    'FontSize', 18, 'Interpreter', 'LaTeX', 'EdgeColor', 'none' );

% Compute RMS error
xiMin = min(bearingVec);
xiMax = max(bearingVec);

tMin = -0.5;
tMax = 0.5;

xiMeasCompare = bearingVec( (bearingVec >= xiMin) & (bearingVec < xiMax) );
fMeasCompare = fVec( (bearingVec >= xiMin) & (bearingVec < xiMax) );

fErrors = 0.*xiMeasCompare; % Initialize

for xiMeasCount = 1 : length( xiMeasCompare )
    
   % Find closest index
   xiLoop_deg = xiMeasCompare( xiMeasCount );
   [~, xiInd] = min( abs( 180.*xi./pi - xiLoop_deg ) );
   
   % Get prediction and estimated frequency at that time.
   fLoop = fMeasCompare( xiMeasCount );
   fPredicted = fw( xiInd );
   
   % Store
   fErrors( xiMeasCount ) = fPredicted - fLoop;
    
end

% Return RMS error
rmsError = rms( fErrors );
