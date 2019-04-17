function [varargout]=grav_stress(varargin)

% grav_stress.m - program to find the stress due to gravity from gravitational
% potential energy (GPE) values.
% 
% Calculates vertically-averaged, gravitationally-induced deviatoric stress.
% 
% Based on Lagrangian multipliers from Flesch et al. (2001) but solving
% stress balance equation in matrix form.
% Uses cartesian coordinates, assumes flat Earth but takes into account
% smaller spacing between lines of longitude at higher latitudes.
% 
% Requires long, lat, and GPE for points in a gridded region.
% xyz-based output file - designed for use with GMT's psvelo - contains:
% Column 1: Longitude
% Column 2: Latitude
% Column 3: Most extensional stress (SHmin)
% Column 4: Most compressional stress (SHmax)
% Column 5: Azimuth of SHmax
% Column 6: Stress magnitude - sqrt(0.5*tau_ij*tau_ij)
% 
% Written by Hamish Hirschberg 2016-17
% 
% Input variables (choose a below combination)
% 1. xyz file or matlab variable containing long, lat, and GPE
% OR
% 1. xyz file or matlab variable containing long, lat, and GPE AND
% 2. name of output file in xyz-based format
% OR
% 1. matlab vector or m x n matrix containing long AND
% 2. matlab vector or m x n matrix containing lat AND
% 3. matlab vector or m x n matrix containing GPE, each with same dimensions
% OR
% 1. - 3. as above AND
% 4. name of output file in xyz-based format
% 
% Output variables (optional)
% 1. m x n x 3 matlab variable containing xx, yy, and xy components of
%    deviatoric gravitational stress

data=[];                    % initialise data
output=[];
rE=6371;                    % radius of the Earth in km

% read in input
if nargin<=2
    % if using one of first two input options
    if ischar(varargin{1})
        data=dlmread(varargin{1});      % read in xyz file
    else
        data=varargin{1};           % read in matlab variable
    end
    if nargin==2
        output=varargin{2};         % output file name
    end
else
    % or if using one of second two input options
    if size(varargin{3},2)==1
        data=[varargin{1},varargin{2},varargin{3}]; % combine data vectors
    else
        % if data input as matrices
        dlo=varargin{1}(1,2)-varargin{1}(1,1);      % longitude spacing (deg)
        dla=varargin{2}(2,1)-varargin{2}(1,1);      % latitude spacing (deg)
        [ny,nx]=size(varargin{3});                  % dimensions of region
        dy=rE*dla*pi/180;                           % colatitude spacing in km
        dx=rE*cosd(varargin{2})*dlo*pi/180;         % longitude spacing in km
        Gamma=varargin{3};                          % read in GPE
    end
    if nargin==4
        output=varargin{4};     % output file name
    end
end

% read in data and define region using max and min lat and long if data
% input as vectors not matrices
if ~isempty(data)
    lomin=min(data(:,1));                       % minimum longitude
    lomax=max(data(:,1));                       % maximum longitude
    dlo=data(2,1)-data(1,1);               % longitude spacing
    nx=abs(round((lomax-lomin)/dlo))+1;              % number of x/longitude points
    
    lamin=min(data(:,2));                        % minimum latitude
    lamax=max(data(:,2));                        % maximum latitude
    dla=data(nx+1,2)-data(1,2);                  % latitude spacing
    ny=abs(round((lamax-lamin)/dla))+1;          % number of y/latitude points
    
    dy=rE*dla*pi/180;                           % colatitude spacing in km
    dx=zeros(ny,nx);                             % azimuthal spacing in km
    Gamma=zeros(ny,nx);                    % matrix of GPE
    
    % create matrix of GPE and dx
    for ii=1:ny
        for j=1:nx
            k=(ii-1)*nx+j;          % index of point
            % account for smaller longitude spacing at higher latitudes
            dx(ii,j)=rE*cosd(data(k,2))*dlo*pi/180;
            Gamma(ii,j)=data(k,3);      % read in data
        end
    end
end

nin=(nx-2)*(ny-2);   % number of points not on the edge
ntot=ny*nx;         % total number of points

% start on actual calculations
% take gradients of Gamma
[ Gamx,Gamy ]=gradient(Gamma,1,dy);
Gamx=Gamx./dx;      % account for spacing

taug=zeros(ny,nx,3);

% find stress field with Lam=(lamth,lamph) at given point
strbal=zeros(2*nin);          % matrix to solve stress balance equation in interior
Gamv=zeros(2*nin,1);            % vector of Gamma derivatives

for ii=2:ny-1
    for j=2:nx-1
        k=(ii-2)*(nx-2)+j-1;

        % coefficients of stress balance equation
        strbal(k,k)=-4/dx(ii,j)^2-1/dy^2;
        strbal(k+nin,k+nin)=-1/dx(ii,j)^2-4/dy^2;
        % gradients of Gamma into vector form
        Gamv(k)=Gamx(ii,j);
        Gamv(k+nin)=Gamy(ii,j);

        % accounting for bottom boundary - k(nb)=j
        if ii~=2
            strbal(k,k-nx+2)=1/(2*dy^2);            % below
            strbal(k+nin,k+nin-nx+2)=2/dy^2;
            if j~=nx-1
                strbal(k,k+nin-nx+3)=-3/(8*dx(ii,j)*dy);      % below right
                strbal(k+nin,k-nx+3)=-3/(8*dx(ii,j)*dy);
            end
            if j~=2
                strbal(k,k+nin-nx+1)=3/(8*dx(ii,j)*dy);       % below left
                strbal(k+nin,k-nx+1)=3/(8*dx(ii,j)*dy);
            end
        end

        % accounting for top boundary - k(nb)=2*nx+ny-1-j
        if ii~=ny-1
            strbal(k,k+nx-2)=1/(2*dy^2);            % above
            strbal(k+nin,k+nin+nx-2)=2/dy^2;
            if j~=nx-1
                strbal(k,k+nin+nx-1)=3/(8*dx(ii,j)*dy);      % above right
                strbal(k+nin,k+nx-1)=3/(8*dx(ii,j)*dy);
            end
            if j~=2
                strbal(k,k+nin+nx-3)=-3/(8*dx(ii,j)*dy);       % above left
                strbal(k+nin,k+nx-3)=-3/(8*dx(ii,j)*dy);
            end
        end

        % accounting for right-hand boundary - k(nb)=nx-1+ii
        if j~=nx-1
            strbal(k,k+1)=2/dx(ii,j)^2;
            strbal(k+nin,k+nin+1)=1/(2*dx(ii,j)^2);
        end

        % accounting for left-hand boundary -
        % k(nb)=2*(nx+ny)-2-ii=nb+2-ii
        if j~=2
            strbal(k,k-1)=2/dx(ii,j)^2;
            strbal(k+nin,k+nin-1)=1/(2*dx(ii,j)^2);
        end
    end
end

% solve equation strbal*lambda=strbc*lambc+Gamv for lambda
lambda=strbal\(Gamv);       % calculation performed here

lamx=zeros(ny,nx);
lamy=zeros(ny,nx);

% converting lambda to matrices of lambda_x and lambda_y,
% in the interior (boundary values all 0)
for ii=2:ny-1
    for j=2:nx-1
        k=(ii-2)*(nx-2)+j-1;
        lamx(ii,j)=lambda(k);
        lamy(ii,j)=lambda(k+nin);
    end
end

% derivatives of lambdas to taus
[ lamxx,lamxy ]=gradient(lamx);
[ lamyx,lamyy ]=gradient(lamy);
taug(:,:,1)=lamxx./dx;              % account for grid spacing
taug(:,:,2)=lamyy/dy;
taug(:,:,3)=(lamxy/dy+lamyx./dx)/2;

% Calculate principal stress axes
paaz=princaxes(taug);
nans=isnan(paaz);       % identify NaN resulting from 0/0
paaz(nans)=0;           % replace NaN with 0

if ~isempty(output)
    % print out principal axes for plotting in GMT
    outmat=zeros(ntot,7);           % matrix to be turned into file
    
    % find square-root of second invariant
    taugii=sqrt(taug(:,:,1).^2+taug(:,:,1).*taug(:,:,2)+...
        taug(:,:,2).^2+taug(:,:,3).^2);
    
    % set long and lat in output matrix
    if isempty(data)
        % if data input as matrices
        for ii=1:ny
            for j=1:nx
                k=(ii-1)*nx+j;
                outmat(k,1)=varargin{1}(ii,j);
                outmat(k,2)=varargin{2}(ii,j);
            end
        end
    else
        % or if data input as vectors
        outmat(:,1:2)=data(:,1:2);
    end
    
    % fill out rest of output matrix
    for ii=1:ny
        for j=1:nx
            k=(ii-1)*nx+j;
            outmat(k,3)=paaz(ii,j,1);   % most extensional stress (SHmin) in MPa
            outmat(k,4)=paaz(ii,j,2);   % most compressive stress (SHmax) in MPa
            outmat(k,5)=paaz(ii,j,3);   % azimuth of SHmax in deg
            outmat(k,6)=taugii(ii,j);   % magnitude of gravitational stress in MPa
        end
    end
    
    dlmwrite(output,outmat,'\t');       % write out output file
end

if nargout>=1
    varargout{1}=taug;          % output stress as m x n x 3 variable
    if nargout>=2
        varargout{2}=Gamma;
    end
end

end



