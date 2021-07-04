function output=lpntrp(N, dly)
%LPNTRP   Function to generate coefficients of the fractional delay interpolator.
%   B = LPNTRP(N, DELAY) designs interpolation filter of the order N > 2 which 
%   delays input signal by DELAY of cycles cycles such as 0 < DELAY < 1. 
%
%   Notice, vector B will have N+1 elements, but overall filter delay will be
%   FlOOR(N/2)+DELAY

% Overall delay
D=floor(N/2)+dly;

output=prod(D-flipud(reshape(repmat(0:N,N,1), N+1,N)).')./prod(repmat(0:N,N+1,1)-repmat((0:N)',1,N+1)+diag(ones(1,N+1)));
    


