function [R1int,R2int,tOut] = InterpolateRotMatrices(R1,R2,t1,t2)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

% convert rotation matrices to euler angles
q1 = rotm2eul(R1);
q2 = rotm2eul(R2);

% output time
tOut = unique([t1 t2]);

% interpolate euler angles
q1_int = interp1(t1,q1,tOut);
q2_int = interp1(t2,q2,tOut);

% convert to rotation matrices again
R1int =  eul2rotm(q1_int);
R2int =  eul2rotm(q2_int);



end

