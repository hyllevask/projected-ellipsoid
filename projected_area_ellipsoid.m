function [Ax,Ay,Az] = projected_area_ellipsoid(Lx,Ly,Lz,phi,theta)
%projected_area_ellipsoid Calculates the projected area in the x,y,z
%directions.
% See https://math.stackexchange.com/questions/2438495/showing-positive-definiteness-in-the-projection-of-ellipsoid
% for explenations.
%------------------------------------------------
%   Inputs:
%   Lx:     Length of unrotated x axis
%   Ly:     Length of unrotated y axis
%   Lz:     Length of unrotated z axis
%   theta:  Azimuthal tilt of the ellipsoid
%   phi:    Polar angle
%
%   Outputs:
%   Ax:     Projected area in x direction
%   Ay:     Projected area in y direction
%   Az: 	Projected area in z direction
%------------------------------------------------

%Create the eigenvectors from the input angles
v1 = [cos(phi)*cos(theta),sin(phi)*cos(theta),-sin(theta)]'; %"rotated x-axis"
v2 = [-sin(phi),cos(phi),0]';    %"rotated y-axis"
v3 = [sin(theta)*cos(phi),sin(theta)*sin(phi),cos(theta)]'; %"rotated z-axis"
V = [v1,v2,v3];

%Create the eigenvalues from the lengths
D = diag([1/Lx^2,1/Ly^2,1/Lz^2]);
%Calculate the eliipsoid on quadratic form.
A = (V*D)\(V);


%Calculate the projected ellipse
Px = Calc_P(A,'x');
Py = Calc_P(A,'y');
Pz = Calc_P(A,'z');

%Calculate the area of the projected ellipses
Ax = Calc_area(Px);
Ay = Calc_area(Py);
Az = Calc_area(Pz);

end
function Area = Calc_area(P)
% The eigenvalues of P containes information about the major/minor axis
[~,PD] = eig(P);
Area = pi/(sqrt(PD(1,1)*PD(2,2)));
end
function P = Calc_P(A,dir)
%Calculates the projected ellipse for the three different directions.
%Calculated from the criteria that dot(nabla(x'Ax - 1),n) = 0,
%where n is the direction of the projection.
if dir == 'x'
    P = [A(2,2)-A(1,2)^2/A(1,1), A(2,3)-A(1,2)*A(1,3)/A(1,1);...
        A(2,3)-A(1,2)*A(1,3)/A(1,1), A(3,3)-A(1,3)^2/A(1,1)];
    
elseif dir == 'y'
    P = [A(1,1)-A(1,2)^2/A(2,2), A(1,3)-A(2,3)*A(1,2)/A(2,2); ...
         A(1,3)-A(2,3)*A(1,2)/A(2,2), A(3,3) - A(2,3)^2/A(3,3)];
    
elseif dir == 'z'
    P = [A(1,1)-A(1,3)^2/A(3,3),A(1,2)-A(1,3)*A(2,3)/A(3,3); ...
         A(1,2)-A(1,3)*A(2,3)/A(3,3), A(2,2)-A(2,3)^2/A(3,3)];
    
else
    error('Invalid Direction')
end

end
