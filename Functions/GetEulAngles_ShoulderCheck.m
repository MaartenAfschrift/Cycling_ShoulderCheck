function [eul] = GetEulAngles_ShoulderCheck(R)
%GetEulAngles_ShoulderCheck Convert Rotation matrix to euler angles for the
%shoulder check task
Rnorm = nan(size(R));
nfr = length(R(1,1,:));
R0 = squeeze(R(:,:,1));
for ifr = 1:nfr
    Rnorm(:,:,ifr) = (R(:,:,ifr)*R0');
end
eul = rotm2eul(Rnorm,'ZYX');
end

