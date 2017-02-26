% lin_comb.m - program to find the linear combination of stress basis
% functions that gives the total stress field that most closely matches the
% orientation of the strain rate field
% Hamish Hirschberg

% taug - nth x nph x 3 array of gravitational stress
% bases - nth x nph x nb x 3 array of nb basis functions
% eps - nth x nph x 3 array of strain rate
% coeffs - nb x 1 vector of coefficients of basis functions
% sins - nth x nph matrix of sines of colatitudes

% % read in modeled velocities from vel_calc.m
% vels=dlmread('vel_GPS_model.xyz','\t');
% u=zeros(ny,nx);
% v=zeros(ny,nx);
% % epsd=zeros(ny,nx,3);
% for ii=1:ny
%     for j=1:nx
%         k=(ii-1)*nx+j;
%         u(ii,j)=vels(k,3)/yearsecs*1e9;          % convert to pm/s=10^-9m/s
%         v(ii,j)=vels(k,4)/yearsecs*1e9;
% %         epsd(ii,j,1)=vels(k,5);
% %         epsd(ii,j,2)=vels(k,6);
% %         epsd(ii,j,3)=vels(k,7);
%     end
% end

% Calculate strain rates from gradients of velocities
epsd=zeros(ny,nx,3);
[ux,uy]=gradient(u);
[vx,vy]=gradient(v);
epsd(:,:,1)=ux./dx;
epsd(:,:,2)=vy./dy;
epsd(:,:,3)=(vx./dx+uy./dy)/2;

% % calculate strain rates through collocation
% [epsxx,epsyy,epsxy]=collocds(los,las,u,v,abs(dy));
% epsd(:,:,1)=epsxx;
% epsd(:,:,2)=epsyy;
% epsd(:,:,3)=epsxy;

% add noise to data to test stability
% epsd=epsd+normrnd(0,0.2,ny,nx,3);       % 0.8 approx 1mm/yr uncertainty in velocity
% taug=taug+normrnd(0,1.5,ny,nx,3);       % 1.5 approx 1 MPa uncertainty in GPE

coeffs=ones(size(bases,3),1);      % create vector of coefficients for basis functions
options=optimset('MaxIter',1e5,'MaxFunEvals',1e6,'TolFun',1e-9,'TolX',1e-8);
func=@(coeffs)totlincomb(taug,bases,epsd,coeffs);        % set up function
[coeffs,out]=fminunc(func,coeffs,options);            % minimise the function

taub=zeros(size(taug));     % set up boundary stresses
for ii=1:length(coeffs)
    taub=taub+squeeze(bases(:,:,ii,:)*coeffs(ii));   % boundary stresses as sum of basis functions
end
taut=taug+taub;      % total stress as sum of gravitational and boundary stresses

% figure('Name','tau_phiphi');
% axis ij;
% pl4=surf(taub(:,:,1));
% figure('Name','tau_thth');
% pl5=surf(taub(:,:,2));
% figure('Name','tau_phith');
% pl6=surf(taub(:,:,3));

% calculate misfit between stress and strain rate
% E, T, etau - indetermediary functions
E=sqrt(2*epsd(:,:,1).^2+2*epsd(:,:,1).*epsd(:,:,2)+2*epsd(:,:,2).^2+2*epsd(:,:,3).^2);
T=sqrt(2*taut(:,:,1).^2+2*taut(:,:,1).*taut(:,:,2)+2*taut(:,:,2).^2+2*taut(:,:,3).^2);
etau=2*epsd(:,:,1).*taut(:,:,1)+taut(:,:,1).*epsd(:,:,2)+taut(:,:,2).*epsd(:,:,1)+...
    2*epsd(:,:,2).*taut(:,:,2)+2*taut(:,:,3).*epsd(:,:,3);

misfit=(1-etau./(E.*T))/2;
tol=1;       % maximum acceptable misfit
% mine=0.2;     % minimum accepted value of E
bad=(isnan(misfit)|misfit>tol|E<mine);     % ignore points that have a misfit greater than tol


% calculate principal axes of boundary stress
[paazb,prinaxb]=princaxes(taub);
prinaxb(repmat(bad,[1,1,4]))=NaN;
paazb(repmat(bad,[1,1,3]))=0;
B=sqrt(2*taub(:,:,1).^2+2*taub(:,:,1).*taub(:,:,2)+2*taub(:,:,2).^2+2*taub(:,:,3).^2);
B(bad)=NaN;

% % plot misfit and principal axes of boundary stresses
% factor=5/(max(prinaxb(:)));
% figure('Name','Misfit and Principal Boundary Stress');
% pl1=contour(misfit);
% ax=gca;
% set(ax,'DataAspectRatio',[abs(dy(1,1)),dx(1,1),dx(1,1)]);
% hold on;
% axis ij;
% quiver(prinaxb(:,:,1)*factor,prinaxb(:,:,2)*factor,0,'Color',[0.5 0.5 0.5]);
% quiver(prinaxb(:,:,3)*factor,prinaxb(:,:,4)*factor,0,'Color',[0.3 0.3 0.3]);
% hold off;

% calculate viscosity
eta=T./E;
eta(bad)=NaN;       % ignore viscosities if misfit too large
T(bad)=NaN;

% calculate principal axes of strain rate
[paaze,prinaxe]=princaxes(epsd);
paaze(repmat(bad,[1,1,3]))=0;      % convert unreliable values to 0

% % plot misfit and principal axes of strain rates
% factor=2/(max(prinaxe(:)));
% figure('Name','Viscosity and Principal Strain Rates');
% pl2=contour(eta,0:2:20);
% ax=gca;
% set(ax,'DataAspectRatio',[abs(dy(1,1)),dx(1,1),dx(1,1)]);
% hold on;
% axis ij;
% quiver(prinaxe(:,:,1)*factor,prinaxe(:,:,2)*factor,0,'Color',[0.5 0.5 0.5]);
% quiver(prinaxe(:,:,3)*factor,prinaxe(:,:,4)*factor,0,'Color',[0.3 0.3 0.3]);
% hold off;

% calculate principal axes of total deviatoric stress
[paazt,pinaxt]=princaxes(taut);
paazt(repmat(bad,[1,1,3]))=0;      % convert NaN to 0

% calculate differential stress, assuming one primary stress vertical
sd=max(paazt(:,:,1),-paazt(:,:,1)-paazt(:,:,2))-min(paazt(:,:,2),-paazt(:,:,1)-paazt(:,:,2));
sd(bad)=NaN;

% print out principal axes for plotting in GMT
outmat=zeros(ntot,7);           % matrix to be turned into file
outmatb=zeros(ntot,8);          % second matrix to be turned into file
outmatc=zeros(ntot,8);          % third matrix to be turned into file

for ii=1:ny
    for j=1:nx
        k=(ii-1)*nx+j;
        l=(ii-1)*((nx-1)*dec^2+dec)+(j-1)*dec+1;
        
        % first file
        outmat(k,1)=data(l,1);
        outmat(k,2)=data(l,2);
        outmat(k,3)=paazb(ii,j,1);
        outmat(k,4)=paazb(ii,j,2);
        outmat(k,5)=paazb(ii,j,3);
        outmat(k,6)=B(ii,j);
        outmat(k,7)=misfit(ii,j);
        
        % second file
        outmatb(k,1)=data(l,1);
        outmatb(k,2)=data(l,2);
        outmatb(k,3)=paaze(ii,j,1);
        outmatb(k,4)=paaze(ii,j,2);
        outmatb(k,5)=paaze(ii,j,3);
        outmatb(k,6)=E(ii,j);
        outmatb(k,7)=u(ii,j)*yearsecs/1e9;
        outmatb(k,8)=v(ii,j)*yearsecs/1e9;
        
        % third file
        outmatc(k,1)=data(l,1);
        outmatc(k,2)=data(l,2);
        outmatc(k,3)=paazt(ii,j,1);
        outmatc(k,4)=paazt(ii,j,2);
        outmatc(k,5)=paazt(ii,j,3);
        outmatc(k,6)=eta(ii,j);
        outmatc(k,7)=T(ii,j);
        outmatc(k,8)=sd(ii,j);
        
    end
end

dlmwrite(strcat('bound_visc_',model,'.txt'),outmat,'\t');
dlmwrite(strcat('strain_mis_',model,'.txt'),outmatb,'\t');
dlmwrite(strcat('total_str_',model,'.txt'),outmatc,'\t');






