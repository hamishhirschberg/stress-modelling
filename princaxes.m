function [paaz,prinax]=princaxes(varargin)
% princaxes - function to calculate principal axes and principal directions
% for an inputted array of components of a symmetric 2x2 tensor, such as
% horizontal stress or strain.
%
% Input (order: xx, yy, xy)
% EITHER
% 1. m x n x 3 array corresponding to the three components for each of
%   m x n points.
% OR
% 1-3. 3 m x n arrays, each corresponding to one of the three components
%
% Output
% 1. m x n x 3 array of principal axes magnitudes (largest and smallest) 
%   and azimuth of smallest principal axis (intended for plotting in GMT)
% 2. m x n x 4 array of principal axes x- and y-components (intended for
%   plotting in Matlab)
%
% Calculations assume x=East, y=North; extension positive.
%
% Written by Hamish Hirschberg

% if inputting components as three separate matrices
if nargin==3
    comps(:,:,1)=varargin{1};
    comps(:,:,2)=varargin{2};
    comps(:,:,3)=varargin{3};
else
    comps=varargin{1};
end

% find dimensions of region
ny=size(comps,1);
nx=size(comps,2);

% Calculate principal stress axes
prinax=zeros(ny,nx,4);          % principal axes for plotting by matlab
paaz=zeros(ny,nx,3);            % principal axes for plotting by gmt
       
% Using analytically determined expression:
% eigval=(comp_xx+comp_yy+/-sqrt((comp_xx-comp_yy)^2+4comp_xy^2))/2
for ii=1:ny
    for j=1:nx
        discr=hypot(comps(ii,j,1)-comps(ii,j,2),2*comps(ii,j,3));
        
        % first eigenvalue
        evla=(comps(ii,j,1)+comps(ii,j,2)+discr)/2;
        paaz(ii,j,1)=evla;
        evca=(comps(ii,j,2)-comps(ii,j,1)+discr)/(2*comps(ii,j,3));
        norma=hypot(1,evca);        % normalise eigenvector
        prinax(ii,j,1)= evla/norma;  % store and weight by eigenvalue
        prinax(ii,j,2)=-evca*evla/norma;
        
        % second eigenvalue
        evla=(comps(ii,j,1)+comps(ii,j,2)-discr)/2;
        paaz(ii,j,2)=evla;
        evca=(comps(ii,j,2)-comps(ii,j,1)-discr)/(2*comps(ii,j,3));
        % azimuth of eigenvector in degrees between 0 and 360
        paaz(ii,j,3)=wrapTo360(atan(1/evca)*180/pi);
        norma=hypot(1,evca);        % normalise eigenvector
        prinax(ii,j,3)= evla/norma;  % store and weight by eigenvalue
        prinax(ii,j,4)=-evca*evla/norma;
    end
end

end