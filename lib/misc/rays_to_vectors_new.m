function [vectors,rayXYZ, ths_aver,ths_max] = ...
						rays_to_vectors_new(rx,ry, geometry)

% (theta-EQUATOR version)
%
%	Calculate 'vectors' matrix for projection geometry in ASTRA Toolbox
%	from lists of ray director (unit-length) components: rx, ry
%	Normalize zenital angles in illumination scenario
%
%
%%%%%% INPUT:
% rx, ry	- coordinates of the illumination directors (unit length!)
%%%%%% <optional>
%
%%%%%% OUTPUT
% vectors   - projection geometry for ASTRA Toolbox
%
%     ( rayX, rayY, rayZ, dX, dY, dZ, uX, uY, uZ, vX, vY, vZ )
%     ray : the ray direction
%     d   : the center of the detector
%     u   : the vector from detector pixel (0,0) to (0,1)
%     v   : the vector from detector pixel (0,0) to (1,0)
%
% rayXYZ    = [rx;ry;rz] (normalized) illumination scenario

rx0 = reshape(rx, 1, []);
ry0 = reshape(ry, 1, []);
if sum(rx0.^2+ry0.^2-1 > 1e-15)>0
	error('Some illumination versors seem to have over-unit length.');
    %disp('Warning: Some illumination versors seem to have over-unit length.');
end
rz0 = real(sqrt(1-rx0.^2-ry0.^2));

rx = rx0;
ry = ry0;
rz = rz0;	
	
ths = asin(sqrt(rx.^2+ry.^2)); % theta-POLE !!!!!!!!!!!!!!!!!!!!!!!
% ths = acos(min(1,sqrt(rx.^2+ry.^2))); % theta-EQUATOR
ths_aver = mean(ths);
ths_max = max(ths);
ths_aver_deg = rad2deg(mean(ths));
ths_max_deg = rad2deg(max(ths));

rayXYZ = [rx; ry; rz];

% rotation matrices
np = length(rx); % number of projections
Rray = zeros(3,3,np); 
firay = atan2(rayXYZ(2,:), rayXYZ(1,:));
%thray = asin(sqrt(rayXYZ(1,:).^2+rayXYZ(2,:).^2)); % theta-POLE
thray = acos(min(1,sqrt(rayXYZ(1,:).^2+rayXYZ(2,:).^2))); % theta-EQUATOR
for p=1:np % R matrices transforming [0;0;1] into rayXYZ(:,s)
%    Rray(:,:,p) = aX1aX2toDCM(firay(p), 0, 'X', false) * ...
%                  aX1aX2toDCM(-thray(p), 0, 'Z', false); % theta-EQUATOR
%%                  aX1aX2toDCM(thray(p), 0, 'Z', false); % theta-POLE
	if 		abs(thray(p)-pi/2)<1e-12; thray(p)=pi/2-1e-12; % avoid singularity
	elseif 	abs(thray(p)+pi/2)<1e-12; thray(p)=-pi/2+1e-12; end
	Rray(:,:,p) = SpinCalc('EA123toDCM',-[0 0 rad2deg(firay(p))],1e-6,false) * ...
			      SpinCalc('EA123toDCM',-[0 rad2deg(-thray(p)) 0],1e-6,false); % theta-EQUATOR
%			      SpinCalc('EA321toDCM',-[0 rad2deg(thray(p)) 0],1e-6,false); % theta-POLE
end

% detector center shifts
dXYZ = zeros(3,np);

% projection plane versors
switch geometry
case 'fixed' % detector always horizontal
    
    uXYZ = repmat([1; 0; 0], [1 np]);
    vXYZ = repmat([0; 1; 0], [1 np]);
    
case 'facing' % detector perpendicular to illumination (follows the rays)
    
    uXYZ = zeros(3,np); % u==x
    vXYZ = zeros(3,np); % v==y
    det_spacing_x = 1;
    det_spacing_y = 1;
    for p=1:np
        % looking from the positive side of Z => illumination angles with + sign
%        uXYZ(:,p) = det_spacing_x * squeeze(Rray(:,:,p)) * [1; 0; 0]; % u: theta-POLE
        uXYZ(:,p) = det_spacing_x * squeeze(Rray(:,:,p)) * [0; 0; -1]; % u: theta-EQUATOR
        vXYZ(:,p) = det_spacing_y * squeeze(Rray(:,:,p)) * [0; 1; 0];
    end
    
otherwise
    error('Wrong option `projections`')
end

% projection vectors
vectors = [rayXYZ' dXYZ' uXYZ' vXYZ'];

end
