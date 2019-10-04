%**************************************************************************
%
% ASA Student Challenge Problem 2019
%
%   Problem 1 - Plots the incident angle of the arrival ray as a function
%   of time, as well as the resulting Doppler shifted frequency as a
%   function of time, for comparison with the provided measurements.
%
%
%        Scott Schoen Jr | Georgia Tech | 30 September 2019
%
%**************************************************************************

clear all
close all
clc

dataFile = '../media/timeVsFreq.txt';

% Load in data
allData = importdata( dataFile );
tMeas = allData( :, 1 );
fMeas = allData(:, 2 );

% Problem parameters
d = 20; % [m]
h = 151; % [m]
v = 123; % [m/s]
f0 = 68; % [Hz]

c1 = 343; % [m/s]
c2 = 1520; % [m/s]

% Define plane position at each point
t = linspace( -1.4, 1.4, 1E3 );
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
        D - (h.*tan(thetaVar) + d.*tan( asin((c2./c1).*sin(thetaVar)) ) );    
    theta1(tCount) = fzero( zeroFunction, theta0 );     
    
end


% Plot theta1 as a function of time
figure();
hold all;

plot( t, 180.*theta1./pi, 'k' );
plot( t, 0.*t + 180.*thetac./pi, '--k', 'LineWidth', 2 );
plot( t, 0.*t - 180.*thetac./pi, '--k', 'LineWidth', 2  );

xlabel( 'Time [s]', 'FontSize', 24 );
ylabel( '$\theta_{1}$ [deg]', 'FontSize', 24 );

grid on;
box off;

% Compute theta2
theta2 = asin( (c2./c1).*sin(theta1) );

% Compute the instantaneous frequency at the surface
fw = f0.*( 1 - (v./c1).*sin( theta1 ) ).^(-1);

% Plot frequency vs time
figure();
hold all;

plot( t, fw, 'k', 'LineWidth', 3 ); 
plot( tMeas, fMeas, 'ko', 'MarkerFaceColor', 'k' );

ylabel( 'Frequency [Hz]', 'FontSize', 24 );
xlabel( 'Time [s]', 'FontSize', 24 );

ylim([60, 80]);
grid on;
box off;

% Compute RMS error
tMin = min(tMeas);
tMax = max(tMeas);

tMeasCompare = tMeas( (tMeas >= tMin) & (tMeas < tMax) );
fMeasCompare = fMeas( (tMeas >= tMin) & (tMeas < tMax) );

fErrors = 0.*tMeasCompare; % Initialize

for tMeasCount = 1 : length( tMeasCompare )
    
   % Find closest index
   tLoop = tMeasCompare( tMeasCount );
   [~, tInd] = min( abs( t - tLoop ) );
   
   % Get prediction and estimated frequency at that time.
   fLoop = fMeasCompare( tMeasCount );
   fPredicted = fw( tInd );
   
   % Store
   fErrors( tMeasCount ) = fPredicted - fLoop;
    
end

