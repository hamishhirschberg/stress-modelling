% grav_mat.m - program to find the stress due to gravity from gravitational
% potential energy (GPE) values.
% Hamish Hirschberg
% Based on Lagrangian multipliers from Flesch et al. (2001) but solving
% stress balance equation in matrix form.
% Uses cartesian coordinates, assumes flat Earth but takes into account
% smaller spacing between lines of longitude at higher latitudes

% read in data - only use location of boundary defined by min and max long
% and lat
data=dlmread('GPE_Vel_gridfile2.xyz',' ');     % read in tab-delimited data
dec=10;         % decimation factor in each dimension

lomin=min(data(:,1));                       % minimum longitude
lomax=max(data(:,1));                       % maximum longitude
dlo=abs(data(1,1)-data(1+dec,1));               % longitude spacing
nx=round((lomax-lomin)/dlo)+1;                     % number of x/longitude points

lamin=min(data(:,2));                        % minimum latitude
lamax=max(data(:,2));                        % maximum latitude
dla=data((nx-1)*dec^2+dec+1,2)-data(1,2);                   % latitude spacing
ny=abs(round((lamax-lamin)/dla))+1;                  % number of y/latitude points

rE=6378;                                    % radius of the Earth in km
yearsecs=31556736;                          % number of seconds in a year
dy=rE*dla*pi/180;                             % colatitude spacing in rad
dx=zeros(ny,nx);                             % azimuthal spacing in rad
nin=(nx-2)*(ny-2);   % number of points not on the edge
ntot=ny*nx;         % total number of points

Gamma=zeros(ny,nx);                     % matrix of gravitational potential energy
u=zeros(ny,nx);     % x-velocity
v=zeros(ny,nx);     % y-velocity

nb=ntot-nin;                        % number of points on boundary

% create matrix of sins of colatitude
for ii=1:ny
    for j=1:nx
        k=(ii-1)*((nx-1)*dec^2+dec)+(j-1)*dec+1;
        % account for smaller longitude spacing at higher latitudes
        dx(ii,j)=rE*cos(data(k,2)*pi/180)*dlo*pi/180;
        Gamma(ii,j)=data(k,3);      % read in data
        u(ii,j)=data(k,4)/yearsecs*1e9;          % x-velocity (East) in pm/s=10^-9m/s
        v(ii,j)=data(k,5)/yearsecs*1e9;          % y-velocity (North) in pm/s
    end
end
% add noise to GPE
% Gamma=Gamma+normrnd(0,5,ny,nx);     % std dev of 1 MPa

% take gradients of Gamma
[ Gamx,Gamy ]=gradient(Gamma,1,dy);
Gamx=Gamx./dx;      % account for spacing

taug=zeros(ny,nx,3);

% find stress field with Lam=(lamth,lamph) at given point
strbal=zeros(2*nin);                % matrix to solve stress balance equation in interior
strbc=zeros(2*nin,2*nb);        % matrix of coefficients of boundary lambda
lambc=zeros(2*nb,1);            % matrix of boundary lambda
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

% solve equation strbal*lambda=strbc*lambc+Gamv for lambda
lambda=strbal\(Gamv);

lamx=zeros(ny,nx);
lamy=zeros(ny,nx);

% converting lambda to matrices of lambda_x and lambda_y,
% in the interior...
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
lamx(2:ny,1)=fliplr(lambc(nb/2+nx:nb));
lamy(2:ny,1)=fliplr(lambc(3*nb/2+nx:2*nb));

% derivatives of lambdas to taus
[ lamxx,lamxy ]=gradient(lamx);
[ lamyx,lamyy ]=gradient(lamy);
taug(:,:,1)=lamxx./dx;              % account for grid spacing
taug(:,:,2)=lamyy/dy;
taug(:,:,3)=(lamxy/dy+lamyx./dx)/2;

deriv=zeros(3*nin,2*nin);           % matrix of derivatives
for ii=2:ny-1
    for j=2:nx-1
        k=(ii-2)*(nx-2)+j-1;
        
        if j~=2
            deriv(k,k-1)=-1/dx(ii,j);
            deriv(k+2*nin,k+nin-1)=-1/(2*dx(ii,j));
        end
        if j~=nx-1
            deriv(k,k+1)=1/dx(ii,j);
            deriv(k+2*nin,k+nin+1)=1/(2*dx(ii,j));
        end
        
        if ii~=2
            deriv(k+nin,k+nin-nx+2)=-1/dy;
            deriv(k+2*nin,k-nx+2)=-1/(2*dy);
        end
        if ii~=ny-1
            deriv(k+nin,k+nin+nx-2)=1/dy;
            deriv(k+2*nin,k+nx-2)=1/(2*dy);
        end
    end
end

% Calculate variances and covariances to estimate uncertainty
a=inv(strbal);      % inverse of matrix of ceoffeicients
cov=a*a';           % covariance matrix
var=diag(cov);      % vector of variances
stdev=zeros(ny,nx);    % matrix of standard deviations, simplify to single value per point
covt=deriv*cov*deriv';  % covariance of stress directly
stdevt=zeros(ny,nx);
for ii=2:ny-1
    for j=2:nx-1
        k=(ii-2)*(nx-2)+j-1;
        stdev(ii,j)=sqrt((var(k)+var(k+nin))/2)/abs(dx(ii,j)*dy);      % account for derivatives
        stdevt(ii,j)=sqrt((covt(k,k)+covt(k+nin,k+nin)+covt(k+2*nin,k+2*nin))/3);
    end
end

% Calculate principal stress axes
prinax=zeros(ny,nx,4);          % principal axes for plotting by matlab
paaz=zeros(ny,nx,3);            % principal axes for plotting by gmt

% Using analytically determined expression:
% eigval=(tauxx+tauyy+/-sqrt((tauxx-tauyy)^2+4tauxy^2))/2
for ii=1:ny
    for j=1:nx
        discr=hypot(taug(ii,j,1)-taug(ii,j,2),2*taug(ii,j,3));
        
        % first eigenvalue
        evla=(taug(ii,j,1)+taug(ii,j,2)+discr)/2;
        paaz(ii,j,1)=evla;
        evca=(taug(ii,j,2)-taug(ii,j,1)+discr)/(2*taug(ii,j,3));
        norma=hypot(1,evca);        % normalise eigenvector
        prinax(ii,j,1)= evla/norma;  % store and weight by eigenvalue
        prinax(ii,j,2)=-evca*evla/norma;
        
        % second eigenvalue
        evla=(taug(ii,j,1)+taug(ii,j,3)-discr)/2;
        paaz(ii,j,2)=evla;
        evca=(taug(ii,j,2)-taug(ii,j,1)-discr)/(2*taug(ii,j,3));
        % azimuth of eigenvector in degrees between 0 and 360
        paaz(ii,j,3)=wrapTo360(atan(1/evca)*180/pi);
        norma=hypot(1,evca);        % normalise eigenvector
        prinax(ii,j,3)= evla/norma;  % store and weight by eigenvalue
        prinax(ii,j,4)=-evca*evla/norma;
    end
end

nans=isnan(paaz);       % identify NaN resulting from 0/0
paaz(nans)=0;           % replace NaN with 0

% plot principal axes
factor=1/(max(prinax(:)));
figure('Name','Gamma and Principal Axes');
pl17=contour(Gamma);
ax=gca;
set(ax,'DataAspectRatio',[abs(dy(1,1)),dx(1,1),dx(1,1)]);
hold on;
axis ij;
quiver(prinax(:,:,1)*factor,prinax(:,:,2)*factor,0,'k');
quiver(prinax(:,:,3)*factor,prinax(:,:,4)*factor,0,'b');
hold off;

% print out principal axes for plotting in GMT
outmat=zeros(ntot,7);           % matrix to be turned into file

% find squareroot of second invariant
taugii=sqrt(2*taug(:,:,1).^2+2*taug(:,:,1).*taug(:,:,2)+2*taug(:,:,2).^2+2*taug(:,:,3).^2);

for ii=1:ny
    for j=1:nx
        k=(ii-1)*nx+j;
        l=(ii-1)*((nx-1)*dec^2+dec)+(j-1)*dec+1;
        outmat(k,1)=data(l,1);
        outmat(k,2)=data(l,2);
        outmat(k,3)=Gamma(ii,j);
        outmat(k,4)=paaz(ii,j,1);
        outmat(k,5)=paaz(ii,j,2);
        outmat(k,6)=paaz(ii,j,3);
        outmat(k,7)=taugii(ii,j);
    end
end

dlmwrite('prin_ax_grav.txt',outmat,'\t');



