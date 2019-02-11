% *************************************************************************
% Driver model setting
% *************************************************************************

% Copyright 2018 The MathWorks, Inc.

% sampling time
Ts = 0.1; % [s]

% Among the parameters set below, whichever are present in Widemann-driver model 
% mask they get overwritten with the mask value.

% vehicle parameters
Ln1 = 4.5; % [m] Legnth, vehicle n-1
vdes = 40; % [m/s] desired speed. 
vmax = 60; % [m/s] maximum speed
% wiedemann model parameters
AXadd = 2.5;
AXmult = 0;
BXadd = 2.5;
BXmult = 0;
EXadd = 1.5;
EXmult = 0;
CXconst = 40;
CXadd = 1.0;
CXmult = 0;
OPDVmult = 0;
OPDVadd = 2.25;
BNULLmult = 0.1;
BMAXmult = 3.5/vmax;
FAKTORVmult = 0;
BMINadd = 20;
BMINmult = 1.5/vmax;
% normally distributed driver dependent parameters
RNDmu = [0.5;0.5;0.5;0.5]; % mean
RNDsigma = [0.0;0.0;0.0;0.0]; % variance
RND1n = RNDmu(1) + randn(1,1)*sqrt(RNDsigma(1));
RND2n = RNDmu(2) + randn(1,1)*sqrt(RNDsigma(2));
RND3n = RNDmu(3) + randn(1,1)*sqrt(RNDsigma(3));
RND4n = RNDmu(4) + randn(1,1)*sqrt(RNDsigma(4));
% normally distributed random number
NRNDmu = 0.5; % mean
NRNDsigma = 0.0; % variance
NRND = NRNDmu(1) + randn(1,1)*sqrt(NRNDsigma(1));

W_SPEC = [Ln1;vdes;vmax;AXadd;AXmult;BXadd;BXmult;EXadd;EXmult;CXconst;CXadd;CXmult;OPDVmult;OPDVadd;BNULLmult;BMAXmult;FAKTORVmult;BMINadd;BMINmult; ...
    RND1n;RND2n;RND3n;RND4n;NRND];