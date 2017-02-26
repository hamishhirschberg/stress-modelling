% basis_fn_mat.m - program to find the basis functions corresponding to
% relative plate motions for a given geometry.
% Hamish Hirschberg
% Based on Lagrangian multipliers from Flesch et al. (2001) but solving
% stress balance equation in matrix form in the interior of the grid and
% setting multipliers equal to given values on boundary
% Uses cartesian coordinates, assumes flat Earth
% Creates a basis function for an entire side

% % read in data - only use location of boundary defined by min and max long
% % and lat
% data=dlmread('../GPEsurface_cr.xyz','\t');       %read tab-delimited input file
% ntot=length(data(:,1));                       % total number of points
% dec=1;
% 
% lomin=min(data(:,1));                       % minimum longitude
% lomax=max(data(:,1));                       % maximum longitude
% dlo=abs(data(1,1)-data(2,1));               % longitude spacing
% nx=(lomax-lomin)/dlo+1;                     % number of x/longitude points
% 
% lamin=min(data(:,2));                        % minimum latitude
% lamax=max(data(:,2));                        % maximum latitude
% ny=ntot/nx;                                 % number of y/latitude points
% dla=(lamax-lamin)/(ny-1);                   % latitude spacing
% 
% dy=50;    %dla*pi/180;                             % colatitude spacing in rad
% dx=50*ones(ny,nx);    %dlo*pi/180;                             % azimuthal spacing in rad
% rE=1;                                    % radius of the Earth in km
% 
% nin=(nx-2)*(ny-2);                  % number of points in the interior of the grid
% nb=ntot-nin;                        % number of points on boundary

azs=data(1:dec:(nx-1)*dec+1,1)*pi/180;                        % vector of azimuths/longitudes
els=data(1:(nx-1)*dec^2+dec:(ny-1)*((nx-1)*dec^2+dec)+1,2)*pi/180;                  % vector of elevations/latitudes

boundx=[1 10 20 30];
% boundy=[1 8 15 23];
% boundx=1:5:nx;
% boundy=[1:5:ny,ny];
boundy=[1 10 20 29];
% boundx=[1 10 11 20 21 30 31 40 41 50 51 59];                   % indices of segment endpoints, top & bottom
% boundy=[1 9 10 18 19 27 28 36 37 45];                   % indices of segment endpoints, left & right
sins=zeros(ny,nx);                        % sin(th)

% create matrix of sins of colatitude
for ii=1:ny
    for j=1:nx
        k=(ii-1)*((nx-1)*dec^2+dec)+(j-1)*dec+1;
        % account for smaller longitude spacing at higher latitudes
        th=(90-data(k,2))*pi/180;
        sins(ii,j)=sin(th);
    end
end

% find stress field with Lam=(lamth,lamph) at given point
strbal=zeros(2*nin);                % matrix to solve stress balance equation in interior
lambc=zeros(2*nb,1);            % matrix of lambda values on the boundary
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
bases=zeros(ny,nx,3*(length(boundy)+length(boundx))+3,3);
m=0;            % index of basis function

omg=eye(3);             % matrix of orthogonal vectors

% loop through omega
for l=1:3
    lama=zeros(2,length(boundx));   % matrix of set lambda values
    % get lambda for segment endpoints on bottom side
    for j=1:length(boundx)
        ii=boundx(j);
        [kx,ky,kz]=sph2cart(azs(ii),els(1),1);     % convert point to cartesian coordinates
        % transformation matrix from cartesian to spherical coordinates, ignoring
        % radius
        transf=[sin(els(1)).*cos(azs(ii)) sin(els(1)).*sin(azs(ii)) -sins(1,ii);
            -sin(azs(ii)) cos(azs(ii)) 0];
        % test with a single omega
        Lamcart=cross(omg(:,l),[kx;ky;kz]);             % take cross product Lam=omg x r
        lama(:,j)=transf*Lamcart;                          % convert Lam to spherical
    end
    
    % apply bottom side
    for c=1:length(boundx)-1
        m=m+1;          % index of this basis function
        a=boundx(c);    % location of endpoints of this segment
        b=boundx(c+1);
    
        lambc=zeros(2*nb,1);
        % interpolate lambdax
        lambc(a:b)=spline(boundx,lama(2,:),a:b);
        % interpolate lambday
        lambc(nb+a:nb+b)=spline(boundx,lama(1,:),a:b);

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
    
    lama=zeros(2,length(boundy));   % matrix of set lambda values
    % get lambda for segment endpoints on right-hand side
    for j=1:length(boundy)
        ii=boundy(j);
        [kx,ky,kz]=sph2cart(azs(nx),els(ii),1);     % convert point to cartesian coordinates
        % transformation matrix from cartesian to spherical coordinates, ignoring
        % radius
        transf=[sin(els(ii)).*cos(azs(nx)) sin(els(ii)).*sin(azs(nx)) -sins(ii,nx);
            -sin(azs(nx)) cos(azs(nx)) 0];
        % test with a single omega
        Lamcart=cross(omg(:,l),[kx;ky;kz]);             % take cross product Lam=omg x r
        lama(:,j)=transf*Lamcart;                          % convert Lam to spherical
    end
    
    % apply right-hand side
    for c=1:length(boundy)-1
        m=m+1;          % index of this basis function
        a=boundy(c);    % location of endpoints of this segment
        b=boundy(c+1);
    
        lambc=zeros(2*nb,1);
        % interpolate lambdax
        lambc(nx-1+a:nx-1+b)=spline(boundy,lama(2,:),a:b);
        % interpolate lambday
        lambc(nb+nx-1+a:nb+nx-1+b)=spline(boundy,lama(1,:),a:b);

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
    
    lama=zeros(2,length(boundx));   % matrix of set lambda values
    % get lambda for segment endpoints on top side
    for j=1:length(boundx)
        ii=boundx(j);
        [kx,ky,kz]=sph2cart(azs(ii),els(ny),1);     % convert point to cartesian coordinates
        % transformation matrix from cartesian to spherical coordinates, ignoring
        % radius
        transf=[sin(els(ny)).*cos(azs(ii)) sin(els(ny)).*sin(azs(ii)) -sins(ny,ii);
            -sin(azs(ii)) cos(azs(ii)) 0];
        % test with a single omega
        Lamcart=cross(omg(:,l),[kx;ky;kz]);             % take cross product Lam=omg x r
        lama(:,j)=transf*Lamcart;                          % convert Lam to spherical
    end
    
    % apply top side
    for c=1:length(boundx)-1
        m=m+1;          % index of this basis function
        a=boundx(c);    % location of endpoints of this segment
        b=boundx(c+1);
    
        lambc=zeros(2*nb,1);
        % interpolate lambdax
        lambc(2*nx+ny-1-b:2*nx+ny-1-a)=fliplr(spline(boundx,lama(2,:),a:b));
        % interpolate lambday
        lambc(nb+2*nx+ny-1-b:nb+2*nx+ny-1-a)=fliplr(spline(boundx,lama(1,:),a:b));

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
    
    lama=zeros(2,length(boundy));   % matrix of set lambda values
    % get lambda for segment endpoints on left-hand side
    for j=1:length(boundy)
        ii=boundy(j);
        [kx,ky,kz]=sph2cart(azs(1),els(ii),1);     % convert point to cartesian coordinates
        % transformation matrix from cartesian to spherical coordinates, ignoring
        % radius
        transf=[sin(els(ii)).*cos(azs(1)) sin(els(ii)).*sin(azs(1)) -sins(ii,1);
            -sin(azs(1)) cos(azs(1)) 0];
        % test with a single omega
        Lamcart=cross(omg(:,l),[kx;ky;kz]);             % take cross product Lam=omg x r
        lama(:,j)=transf*Lamcart;                          % convert Lam to spherical
    end
    
    % apply left-hand side
    for c=1:length(boundy)-1
        m=m+1;          % index of this basis function
        a=boundy(c);    % location of endpoints of this segment
        b=boundy(c+1);
    
        lambc=zeros(2*nb,1);
        % treat point (1,1) differently
        if a~=1
            % interpolate lambdax
            lambc(nb+2-b:nb+2-a)=fliplr(spline(boundy,lama(2,:),a:b));
            % interpolate lambday
            lambc(2*nb+2-b:2*nb+2-a)=fliplr(spline(boundy,lama(1,:),a:b));
        else
            % interpolate lambdax
            lambc(nb+2-b:nb)=fliplr(spline(boundy,lama(2,:),a+1:b));
            lambc(1)=spline(boundy,lama(2,:),1);
            % interpolate lambday
            lambc(2*nb+2-b:2*nb)=fliplr(spline(boundy,lama(1,:),a+1:b));
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

% figure('Name','tau_xx');
% pl4=surf(bases(:,:,1,1));
% figure('Name','tau_yy');
% pl5=surf(bases(:,:,1,2));
% figure('Name','tau_xy');
% pl6=surf(bases(:,:,1,3));






