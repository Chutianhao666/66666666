function shift = dopplerShiftCircularOrbit(el,hs,hg,freq,time)
%dopplerShiftCircularOrbit Doppler shift due to satellite movement in circular orbit
%   SHIFT = dopplerShiftCircularOrbit(EL,HS,HG,FREQ) returns the Doppler
%   shift in Hz due to the satellite moving in a circular orbit with the
%   specified elevation angle EL, satellite altitude HS, ground station
%   altitude HG, and satellite carrier frequency FREQ. When the elevation
%   angle is a vector of length NumEL, the output SHIFT is a column vector
%   of the same length. Each row of SHIFT represents the Doppler shift for
%   the corresponding satellite elevation angle. HS, HG, and FREQ are real
%   scalars. The elevation angle is in degrees and the altitudes of
%   satellite and ground station are in meters. This function assumes Earth
%   is spherical, ground station is static, and ignores Earth rotation
%   rate.
%
%   SHIFT = dopplerShiftCircularOrbit(EL,HS,HG,FREQ,TIME) returns the
%   Doppler shift in Hz due to satellite moving in circular orbit at
%   multiple time instances TIME. In this case, the inputs EL, HS, and HG,
%   are treated as the initial values at time equal to 0 second. When TIME
%   is a vector of length NumTIME, SHIFT is a matrix of dimensions
%   NumEL-by-NumTIME. This function assumes that the ground station is
%   located at the North Pole (positive Z-axis) and the satellite starts
%   from the initial input elevation angle in YZ-plane. The satellite is
%   moving in clockwise direction in the circular orbit.
%
%   % Example 1:
%   % Calculate the Doppler shift of a satellite moving in a circular orbit
%   % at an altitude of 600 km and having an elevation angle of 50 degrees
%   % with the ground station on Earth. Assume the satellite carrier
%   % frequency is 2 GHz.
%
%   el = 50;    % degrees
%   hs = 600e3; % meters
%   hg = 0;     % meters
%   freq = 2e9; % Hz
%
%   % Calculate the Doppler shift
%   shift = dopplerShiftCircularOrbit(el,hs,hg,freq)
%
%   % Example 2:
%   % Plot the Doppler shift as a function of elevation angle for a
%   % satellite moving in a circular orbit at an altitude of 10000 km.
%   % Assume the ground station altitude is 120 m and the satellite carrier
%   % frequency is 20 GHz.
%
%   el = 0:90;    % degrees
%   hs = 10000e3; % meters
%   hg = 120;     % meters
%   freq = 20e9;  % Hz
%
%   % Calculate the Doppler shift
%   shift = dopplerShiftCircularOrbit(el,hs,hg,freq);
%
%   % Plot the Doppler shift as function of elevation angle
%   figure
%   plot(el,shift,"-*")
%   title("Doppler Shift vs Elevation Angle")
%   xlabel("Elevation Angle (degrees)")
%   ylabel("Doppler Shift (Hz)")
%   grid on
%
%   % Example 3:
%   % Visualize the variation of Doppler shift for one orbital period
%   % of satellite with the initial elevation angle as 45 degrees. Assume
%   % the satellite altitude is 1500 km and the satellite carrier frequency
%   % is 5 GHz.
%
%   el = 45;       % degrees
%   hs = 1500e3;   % meters
%   hg = 0;        % meters
%   freq = 5e9;    % Hz
%
%   % For the specified satellite altitude of 1500 km, the orbital time
%   % period is 6949.518 seconds. To cover one orbital time period, set the
%   % maximum time instance to 6950 seconds.
%   time = 0:6950; % seconds
%
%   % Calculate Doppler shift for the specified time instances
%   shift = dopplerShiftCircularOrbit(el,hs,hg,freq,time);
%
%   % Plot the Doppler shift as a function of time
%   figure
%   plot(time,shift)
%   title("Doppler Shift vs Time")
%   xlabel("Time (seconds)")
%   ylabel("Doppler Shift (Hz)")
%   grid on
%
%   See also slantRangeCircularOrbit.

%   Copyright 2023 The MathWorks, Inc.

%#codegen

    arguments
        el (1,:) {mustBeFloat, mustBeReal, mustBeFinite, mustBeNonempty}
        hs (1,1) {mustBeFloat, mustBeReal, mustBeFinite, mustBeNonempty, mustBePositive}
        hg (1,1) {mustBeFloat, mustBeReal, mustBeFinite, mustBeNonempty, mustBeNonnegative, mustBeLessThan(hg,hs)}
        freq (1,1) {mustBeFloat, mustBeReal, mustBeFinite, mustBeNonempty}
        time (1,:) {mustBeFloat, mustBeReal, mustBeFinite, mustBeNonempty} = 0
    end

    % Get the statistics of circular orbit
    info = satcom.internal.circularOrbitStats(el,hs,hg,time,freq);

    % Get the Doppler shift
    shift = info.DopplerShift;

end