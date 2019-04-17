function [bases]=basis_fns(varargin)

% basis_fns.m - program to find the stress basis functions corresponding to
% relative plate motions for a given geometry.
%
% Written by Hamish Hirschberg 2016-17
% 
% Based on Lagrangian multipliers from Flesch et al. (2001) but solving
% stress balance equation in matrix form in the interior of the grid and
% setting multipliers equal to given values on boundary.
% Uses cartesian coordinates, assumes flat Earth but takes into account
% smaller spacing of lines of longitude at higher latitudes.
% Creates a basis function for each space between boundary grid points
% 
% Input: EITHER
% 1. region for basis functions as xy(z) file or matlab variable
% OR
% 1. matlab vector or m x n matrix containing long AND
% 2. matlab vector or m x n matrix containing lat of region
% 
% Output
% 1. m x n x k x 3 matlab variable containing k=6*(m+n)-9 stress basis
%    functions

data=[];        % initialise data
rE=6378;        % radius of the Earth in km

if nargin == 1
    % if inputting long and lat together
    if ischar(varargin{1})
        data=dlmread(varargin{1});      % if reading from xy(z) file
    else
        data=varargin{1};           % if using matlab variable
    end
else
    % if inputting long and lat separately
    if size(varargin{1},2)==1
        data=[varargin{1},varargin{2}];     % if using vectors
    else
        % if data input as matrices
        dlo=varargin{1}(1,2)-varargin{1}(1,1);      % longitude spacing (deg)
        dla=varargin{2}(2,1)-varargin{2}(1,1);      % latitude spacing (deg)
        [ny,nx]=size(varargin{3});                  % dimensions of region
        dy=rE*dla*pi/180;                           % colatitude spacing in km
        dx=rE*cosd(varargin{2})*dlo*pi/180;         % longitude spacing in km
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
            k=(ii-1)*(nx-1)+j;          % index of point
            % account for smaller longitude spacing at higher latitudes
            dx(ii,j)=rE*cosd(data(k,2))*dlo*pi/180;
            Gamma(ii,j)=data(k,3);      % read in data
        end
    end
end

nin=(nx-2)*(ny-2);   % number of points not on the edge
ntot=ny*nx;         % total number of points
nb=ntot-nin;        % number of points on boundary

azs=data(1:nx,1)*pi/180;                        % vector of azimuths/longitudes
els=data(1:nx:ntot,2)*pi/180;                  % vector of elevations/latitudes

boundx=1:nx;
boundy=1:ny;

% find stress field with Lam=(lamth,lamph) at given point
strbal=zeros(2*nin);                % matrix to solve stress balance equation in interior
strbc=zeros(2*nin,2*nb);        % matrix of coefficients of bounary lambda

% set up matrices of coefficients representing 2nd derivatives of lambda
for ii=2:ny-1
    for j=2:nx-1
        k=(ii-2)*(nx-2)+j-1;

        % coefficients of stress balance equation
        strbal(k,k)=-4/dx(ii,j)^2-1/dy^2;
        strbal(k+nin,k+nin)=-1/dx(ii,j)^2-4/dy^2;

        % accounting for bottom boundary - k(nb)=j
        if ii==2
            strbc(k,j)=-1/(2*dy^2);
            strbc(k,nb+j+1)=3/(8*dx(ii,j)*dy);            % below right
            strbc(k,nb+j-1)=-3/(8*dx(ii,j)*dy);           % below left
            strbc(k+nin,j+1)=3/(8*dx(ii,j)*dy);
            strbc(k+nin,j-1)=-3/(8*dx(ii,j)*dy);
            strbc(k+nin,j+nb)=-2/dy^2;
        else
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
        if ii==ny-1
            strbc(k,2*nx+ny-1-j)=-1/(2*dy^2);
            strbc(k,nb+2*nx+ny-j)=3/(8*dx(ii,j)*dy);      % above left
            strbc(k,nb+2*nx+ny-2-j)=-3/(8*dx(ii,j)*dy);   % above right
            strbc(k+nin,2*nx+ny-j)=3/(8*dx(ii,j)*dy);
            strbc(k+nin,2*nx+ny-2-j)=-3/(8*dx(ii,j)*dy);
            strbc(k+nin,nb+2*nx+ny-1-j)=-2/dy^2;
        else
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
        if j==nx-1
            strbc(k,nx-1+ii)=-2/dx(ii,j)^2;
            strbc(k,nb+nx+ii)=-3/(8*dx(ii,j)*dy);         % above right
            strbc(k,nb+nx-2+ii)=3/(8*dx(ii,j)*dy);        % below right
            strbc(k+nin,nx+ii)=-3/(8*dx(ii,j)*dy);
            strbc(k+nin,nx-2+ii)=3/(8*dx(ii,j)*dy);
            strbc(k+nin,nb+nx-1+ii)=-1/(2*dx(ii,j)^2);
        else
            strbal(k,k+1)=2/dx(ii,j)^2;
            strbal(k+nin,k+nin+1)=1/(2*dx(ii,j)^2);
        end

        % accounting for left-hand boundary -
        % k(nb)=2*(nx+ny)-2-ii=nb+2-ii
        if j==2
            strbc(k,nb+2-ii)=-2/dx(ii,j)^2;
            strbc(k,2*nb+1-ii)=3/(8*dx(ii,j)*dy);     % above left
            strbc(k+nin,nb+1-ii)=3/(8*dx(ii,j)*dy);
            if ii~=2
                strbc(k,2*nb+3-ii)=-3/(8*dx(ii,j)*dy);    % below left
                strbc(k+nin,nb+3-ii)=-3/(8*dx(ii,j)*dy);
            end
            strbc(k+nin,2*nb+2-ii)=-1/(2*dx(ii,j)^2);
        else
            strbal(k,k-1)=2/dx(ii,j)^2;
            strbal(k+nin,k+nin-1)=1/(2*dx(ii,j)^2);
        end
    end
end

% initialise arrays of stress
bases=zeros(ny,nx,6*(nx+ny)-9,3);
m=0;            % index of basis function

omg=eye(3);             % matrix of orthogonal vectors

% loop through omega
for l=1:3
    lama=zeros(2,nx);   % matrix of set lambda values
    % get lambda for segment endpoints on bottom side
    for ii=1:nx
        [kx,ky,kz]=sph2cart(azs(ii),els(1),1);     % convert point to cartesian coordinates
        % transformation matrix from cartesian to spherical coordinates, ignoring
        % radius
        transf=[sin(els(1)).*cos(azs(ii)) sin(els(1)).*sin(azs(ii)) -cos(els(1));
            -sin(azs(ii)) cos(azs(ii)) 0];
        % test with a single omega
        Lamcart=cross(omg(:,l),[kx;ky;kz]);             % take cross product Lam=omg x r
        lama(:,ii)=transf*Lamcart;                          % convert Lam to spherical
    end
    
    % apply bottom side
    for c=1:nx-1
        m=m+1;          % index of this basis function
    
        lambc=zeros(2*nb,1);
        % interpolate lambdax
        lambc(c:c+1)=spline(boundx,lama(2,:),c:c+1);
        % interpolate lambday
        lambc(nb+c:nb+c+1)=spline(boundx,lama(1,:),c:c+1);

        % solve equation strbal*lambda=strbc*lambc for lambda
        lambda=strbal\(strbc*lambc);

        lamx=zeros(ny,nx);
        lamy=zeros(ny,nx);

        % converting lambda to matrices of lambda_x and lambda_y,
        % in the interior,
        for ii=2:ny-1
            for j=2:nx-1
                k=(ii-2)*(nx-2)+j-1;
                lamx(ii,j)=lambda(k);
                lamy(ii,j)=lambda(k+nin);
            end
        end

        % and on the boundary.
        lamx(1,1:nx-1)=lambc(1:nx-1);
        lamy(1,1:nx-1)=lambc(nb+1:nb+nx-1);
        lamx(1:ny-1,nx)=lambc(nx:nb/2);
        lamy(1:ny-1,nx)=lambc(nb+nx:3*nb/2);
        lamx(ny,2:nx)=flipud(lambc(nb/2+1:nb/2+nx-1));
        lamy(ny,2:nx)=flipud(lambc(3*nb/2+1:3*nb/2+nx-1));
        lamx(2:ny,1)=flipud(lambc(nb/2+nx:nb));
        lamy(2:ny,1)=flipud(lambc(3*nb/2+nx:2*nb));

        % derivatives of lambdas to taus
        [ lamxx,lamxy ]=gradient(lamx);
        [ lamyx,lamyy ]=gradient(lamy);
        bases(:,:,m,1)=lamxx./dx;              % account for grid spacing
        bases(:,:,m,2)=lamyy/dy;
        bases(:,:,m,3)=(lamxy/dy+lamyx./dx)/2;
    end
    
    lama=zeros(2,ny);   % matrix of set lambda values
    % get lambda for segment endpoints on right-hand side
    for ii=1:ny
        [kx,ky,kz]=sph2cart(azs(nx),els(ii),1);     % convert point to cartesian coordinates
        % transformation matrix from cartesian to spherical coordinates, ignoring
        % radius
        transf=[sin(els(ii)).*cos(azs(nx)) sin(els(ii)).*sin(azs(nx)) -cos(els(ii));
            -sin(azs(nx)) cos(azs(nx)) 0];
        % test with a single omega
        Lamcart=cross(omg(:,l),[kx;ky;kz]);             % take cross product Lam=omg x r
        lama(:,ii)=transf*Lamcart;                          % convert Lam to spherical
    end
    
    % apply right-hand side
    for c=1:ny-1
        m=m+1;          % index of this basis function
    
        lambc=zeros(2*nb,1);
        % interpolate lambdax
        lambc(nx-1+c:nx+c)=spline(boundy,lama(2,:),c:c+1);
        % interpolate lambday
        lambc(nb+nx-1+c:nb+nx+c)=spline(boundy,lama(1,:),c:c+1);

        % solve equation strbal*lambda=strbc*lambc for lambda
        lambda=strbal\(strbc*lambc);

        lamx=zeros(ny,nx);
        lamy=zeros(ny,nx);

        % converting lambda to matrices of lambda_x and lambda_y,
        % in the interior,
        for ii=2:ny-1
            for j=2:nx-1
                k=(ii-2)*(nx-2)+j-1;
                lamx(ii,j)=lambda(k);
                lamy(ii,j)=lambda(k+nin);
            end
        end

        % and on the boundary.
        lamx(1,1:nx-1)=lambc(1:nx-1);
        lamy(1,1:nx-1)=lambc(nb+1:nb+nx-1);
        lamx(1:ny-1,nx)=lambc(nx:nb/2);
        lamy(1:ny-1,nx)=lambc(nb+nx:3*nb/2);
        lamx(ny,2:nx)=flipud(lambc(nb/2+1:nb/2+nx-1));
        lamy(ny,2:nx)=flipud(lambc(3*nb/2+1:3*nb/2+nx-1));
        lamx(2:ny,1)=flipud(lambc(nb/2+nx:nb));
        lamy(2:ny,1)=flipud(lambc(3*nb/2+nx:2*nb));

        % derivatives of lambdas to taus
        [ lamxx,lamxy ]=gradient(lamx);
        [ lamyx,lamyy ]=gradient(lamy);
        bases(:,:,m,1)=lamxx./dx;              % account for grid spacing
        bases(:,:,m,2)=lamyy/dy;
        bases(:,:,m,3)=(lamxy/dy+lamyx./dx)/2;
    end
    
    lama=zeros(2,nx);   % matrix of set lambda values
    % get lambda for segment endpoints on top side
    for ii=1:nx
        [kx,ky,kz]=sph2cart(azs(ii),els(ny),1);     % convert point to cartesian coordinates
        % transformation matrix from cartesian to spherical coordinates, ignoring
        % radius
        transf=[sin(els(ny)).*cos(azs(ii)) sin(els(ny)).*sin(azs(ii)) -cos(els(ny));
            -sin(azs(ii)) cos(azs(ii)) 0];
        % test with a single omega
        Lamcart=cross(omg(:,l),[kx;ky;kz]);             % take cross product Lam=omg x r
        lama(:,j)=transf*Lamcart;                          % convert Lam to spherical
    end
    
    % apply top side
    for c=1:nx-1
        m=m+1;          % index of this basis function
    
        lambc=zeros(2*nb,1);
        % interpolate lambdax
        lambc(2*nx+ny-2-c:2*nx+ny-1-c)=fliplr(spline(boundx,lama(2,:),c:c+1));
        % interpolate lambday
        lambc(nb+2*nx+ny-2-c:nb+2*nx+ny-1-c)=fliplr(spline(boundx,lama(1,:),c:c+1));

        % solve equation strbal*lambda=strbc*lambc for lambda
        lambda=strbal\(strbc*lambc);

        lamx=zeros(ny,nx);
        lamy=zeros(ny,nx);

        % converting lambda to matrices of lambda_x and lambda_y,
        % in the interior,
        for ii=2:ny-1
            for j=2:nx-1
                k=(ii-2)*(nx-2)+j-1;
                lamx(ii,j)=lambda(k);
                lamy(ii,j)=lambda(k+nin);
            end
        end

        % and on the boundary.
        lamx(1,1:nx-1)=lambc(1:nx-1);
        lamy(1,1:nx-1)=lambc(nb+1:nb+nx-1);
        lamx(1:ny-1,nx)=lambc(nx:nb/2);
        lamy(1:ny-1,nx)=lambc(nb+nx:3*nb/2);
        lamx(ny,2:nx)=flipud(lambc(nb/2+1:nb/2+nx-1));
        lamy(ny,2:nx)=flipud(lambc(3*nb/2+1:3*nb/2+nx-1));
        lamx(2:ny,1)=flipud(lambc(nb/2+nx:nb));
        lamy(2:ny,1)=flipud(lambc(3*nb/2+nx:2*nb));

        % derivatives of lambdas to taus
        [ lamxx,lamxy ]=gradient(lamx);
        [ lamyx,lamyy ]=gradient(lamy);
        bases(:,:,m,1)=lamxx./dx;              % account for grid spacing
        bases(:,:,m,2)=lamyy/dy;
        bases(:,:,m,3)=(lamxy/dy+lamyx./dx)/2;
    end
    
    lama=zeros(2,ny);   % matrix of set lambda values
    % get lambda for segment endpoints on left-hand side
    for ii=1:ny
        [kx,ky,kz]=sph2cart(azs(1),els(ii),1);     % convert point to cartesian coordinates
        % transformation matrix from cartesian to spherical coordinates, ignoring
        % radius
        transf=[sin(els(ii)).*cos(azs(1)) sin(els(ii)).*sin(azs(1)) -cos(els(ii));
            -sin(azs(1)) cos(azs(1)) 0];
        % test with a single omega
        Lamcart=cross(omg(:,l),[kx;ky;kz]);             % take cross product Lam=omg x r
        lama(:,ii)=transf*Lamcart;                          % convert Lam to spherical
    end
    
    % apply left-hand side
    for c=1:ny-1
        m=m+1;          % index of this basis function
    
        lambc=zeros(2*nb,1);
        % treat point (1,1) differently
        if c~=1
            % interpolate lambdax
            lambc(nb+1-c:nb+2-c)=fliplr(spline(boundy,lama(2,:),c:c+1));
            % interpolate lambday
            lambc(2*nb+1-c:2*nb+2-c)=fliplr(spline(boundy,lama(1,:),c:c+1));
        else
            % interpolate lambdax
            lambc(nb)=fliplr(spline(boundy,lama(2,:),c));
            lambc(1)=spline(boundy,lama(2,:),1);
            % interpolate lambday
            lambc(2*nb)=fliplr(spline(boundy,lama(1,:),c));
            lambc(nb+1)=spline(boundy,lama(1,:),1);
        end

        % solve equation strbal*lambda=strbc*lambc for lambda
        lambda=strbal\(strbc*lambc);

        lamx=zeros(ny,nx);
        lamy=zeros(ny,nx);

        % converting lambda to matrices of lambda_x and lambda_y,
        % in the interior,
        for ii=2:ny-1
            for j=2:nx-1
                k=(ii-2)*(nx-2)+j-1;
                lamx(ii,j)=lambda(k);
                lamy(ii,j)=lambda(k+nin);
            end
        end

        % and on the boundary.
        lamx(1,1:nx-1)=lambc(1:nx-1);
        lamy(1,1:nx-1)=lambc(nb+1:nb+nx-1);
        lamx(1:ny-1,nx)=lambc(nx:nb/2);
        lamy(1:ny-1,nx)=lambc(nb+nx:3*nb/2);
        lamx(ny,2:nx)=flipud(lambc(nb/2+1:nb/2+nx-1));
        lamy(ny,2:nx)=flipud(lambc(3*nb/2+1:3*nb/2+nx-1));
        lamx(2:ny,1)=flipud(lambc(nb/2+nx:nb));
        lamy(2:ny,1)=flipud(lambc(3*nb/2+nx:2*nb));

        % derivatives of lambdas to taus
        [ lamxx,lamxy ]=gradient(lamx);
        [ lamyx,lamyy ]=gradient(lamy);
        bases(:,:,m,1)=lamxx./dx;              % account for grid spacing
        bases(:,:,m,2)=lamyy/dy;
        bases(:,:,m,3)=(lamxy/dy+lamyx./dx)/2;
    end
end

% allow for constant values of tau - one basis each
m=m+1;
bases(:,:,m,1)=ones(ny,nx);
m=m+1;
bases(:,:,m,2)=ones(ny,nx);
m=m+1;
bases(:,:,m,3)=ones(ny,nx);

end