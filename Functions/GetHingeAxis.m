function [R,n_steer,n_frame,qAxis] = GetHingeAxis(R_frame,R_steer,t,dt,varargin)
%GetHingeAxis Computes the rotation axis of a hinge joint
%   Input arguments:
%       (1) R_frame: rotation matrix of body 1 (frame) in world
%       (2) R_steer: rotation matri of body 2 (steer) in world
%       (3) t: time vector
%       (4) dt: time span of callibration movement
%       (5) optional inputs:
%           (a) BoolPlot figures (default is 0)
%   Output arguments:
%       (1) R: rotation matrix to compute rotation of frame w.r.t. steer
%       align the hinge (x euler angle)
%       (2) n_steer: vector with axis in the steer frame
%       (3) n_frame: vector with axis in the bike frame
%       (4) qAxis: euler angles of rotation around functional axis (only
%           when Boolplot is true)



BoolPlot = false;
if ~isempty(varargin)
    BoolPlot = varargin{1};
end

% select time 
iSel = t>=dt(1) & t<=dt(2);
tsel = t(iSel);
nfr = length(tsel);

% adapt input
for i =1:nfr
    R_frame(:,:,i) =R_frame(:,:,i)*R_frame(:,:,1)'; %
    R_steer(:,:,i) =R_steer(:,:,i)*R_steer(:,:,1)'; %
end


% get the data in the time span
R_fr = R_frame(:,:,iSel);
R_st = R_steer(:,:,iSel);


% get rotation matrix with "the steer" expressed in "bike frame" coordinate
% system
Rst_fr = nan(3,3,nfr);
for i =1:nfr
    Rst_fr(:,:,i) = R_fr(:,:,i)'*R_st(:,:,i); %
end
% Rst_fr = Rst_fr*Rst_fr(:,:,1)';

% Create matrix A as in Zehr2007 (equation 2) with both the steer and the
% bike frame expressed in the world coordinate system
A = zeros(nfr*3,6);
for i = 1:nfr
    A((i-1)*3+1 : i*3 , 1:3) = R_st(:,:,i); % rotation matrix steer
    A((i-1)*3+1 : i*3 , 4:6) = -R_fr(:,:,i); % rotation matrix  frame
end

% singular value decomposition on matrix A
[~, ~, V]=svd(A,0); % singular value decomposition (A = U*S*V')

% get the orientation in steer and frame
n_steer = V(1:3,6)/norm(V(1:3,6)); %orientation in steer
n_frame = V(4:6,6)/norm(V(4:6,6)); %orientation in frame

%Compute the residuals
errorFrame = nan(nfr,3);
for i=1:nfr
    errorFrame(i,:)  = eye(3)*n_frame -  Rst_fr(:,:,i)*n_steer;
end
if BoolPlot
    figure();
    set(gcf,'Position',[341   410   929   387]);
    subplot(2,2,1)
    plot(tsel,errorFrame);
    ylabel('Errors [radians ?]');
    xlabel('Time [s]');
end

% adapt coordinate system bike frame such that x-axis aligns with the
% rotation axis
z = cross(n_frame,[0 1 0]);
y = cross(z,n_frame);
R =  [n_frame'; y; z];

%% compute steer angle when bool plot is true
if BoolPlot        
    % compute steer angle
    R_Axis = nan(3,3,nfr);  
    R0 = (R*Rst_fr(:,:,1))';
    for i =1:nfr
        R_Axis(:,:,i) = (R*Rst_fr(:,:,i))*R0; % axis w.r.t. initial frame
    end
    qAxis = rotm2eul(R_Axis,'XYZ');
    subplot(2,2,2);
    plot(tsel,qAxis);
    legend('x-steer','y','z');
    xlabel('Time [s]');
end
if BoolPlot
    % get rotation matrix with "the steer" expressed in "bike frame" coordinate
    % system
    nfr = length(t);
    Rst_fr = nan(3,3,nfr);
    for i =1:nfr
        Rst_fr(:,:,i) = R_frame(:,:,i)'*R_steer(:,:,i); %
    end
    
    % get average orientation of steer w.r.t. frame
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
    subplot(2,2,3:4);
    plot(t,qAxis);
    legend('x-steer','y','z');
    xlabel('time [s]');
    ylabel('Angle [rad]');
    title('full datafile');
else
    qAxis = [];
end








end

