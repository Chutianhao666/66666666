function [fHz, rangeM] = railwayDoppler(motion, carrierHz, lightSpeed)
% x is signed along-track position from closest approach to the base station.
% v is signed dx/dt. Negative velocity at x>0 means approaching, positive CFO.
assert(motion.speedKmh>=0 && any(motion.direction==[-1 1]));
v = motion.direction*motion.speedKmh/3.6;
x = motion.alongTrackM;
rangeM = hypot(x,motion.trackOffsetM);
assert(rangeM>0 && motion.trackOffsetM>=0);
fHz = -carrierHz/lightSpeed*v*x/rangeM;
end
