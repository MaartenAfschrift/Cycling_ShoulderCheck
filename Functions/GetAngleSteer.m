function [qAxis] = GetAngleSteer(R_frame,R_steer,R)
%GetAngleSteer Computes the rotation of body 1 (R_frame) w.r.t. R_steer
%(body 2) in the coordinate system R (with x-axis the hinge joint
%identified in the callibration movement using GetHingeAxis).
%   Input arguments:
%       (1) R_frame: rotation matrix of body 1 (frame) in world
%       (2) R_steer: rotation matri of body 2 (steer) in world
%       (3) R: rotation matrix to compute rotation of frame w.r.t. steer
%       align the hinge (x euler angle)
%   Output arguments:
%       (1) euler angles of rotation R_frame w.r.t. R_steer in the coordinate
%           system R (with x-axis the hinge joint).

% get rotation matrix with "the steer" expressed in "bike frame" coordinate
% system
[~,~,nfr] = size(R_frame);
Rst_fr = nan(3,3,nfr);
for i =1:nfr
    Rst_fr(:,:,i) = R_frame(:,:,i)'*R_steer(:,:,i); %
end

% compute average orientation steer in bike frame
qSteer = rotm2eul(Rst_fr);
qMean = nanmean(qSteer);
RMean =  eul2rotm(qMean);

% in new coordinate system
R_Axis = nan(3,3,nfr);
R0 = (R*RMean)';
for i =1:nfr
    R_Axis(:,:,i) = (R*Rst_fr(:,:,i))*R0;
end

% compute an plot euler angles
qAxis = rotm2eul(R_Axis,'XYZ');


end

