% stress_model.m - program to calculate the vertically averaged deviatoric
% stress following the method of Flesch et al. (2001). This program does
% not do the work itself, consolidating other programs which it calls.
% Hamish Hirschberg

clear
% Using fault motions or gps data?
model='flt';          % calculate the fault motion model
% model='gps';          % calculate the gps model

% % if using GPS model
% if strcmp(model,'gps')
%     vel_model          % run GPS velocity model
% end

% run calculation of gravitational stress
grav_mat_in

% run calculation of basis functions
basis_fn_mat_many

% if using GPS
if strcmp(model,'gps')
    % read in modelled velocities made using GMT surface
    vels=dlmread('NZ_vel_model.xyz');
    u=zeros(ny,nx);
    v=zeros(ny,nx);
    for ii=1:ny
        for j=1:nx
            k=(ii-1)*nx+j;
            u(ii,j)=vels(k,3)/yearsecs*1e9;          % convert to pm/s=10^-9m/s
            v(ii,j)=vels(k,4)/yearsecs*1e9;
        end
    end
end
mine=0.2;

% calculate boundary stresses and put it all together
lin_comb

regs=dlmread('regions.csv',',');        % read in regions matrix
tab=zeros(5,11);             % summary matrix
shmax=squeeze(paazt(:,:,3));        % just the max compressional stress
for ii=1:6
    A=(regs==ii&~bad);          % points in region
    tab(ii,1)=ii;         % region number
    tab(ii,2)=sum(A(:));      % number of points in region
    tab(ii,3)=mean(E(A(:)))*yearsecs/1e8;      % mean strain rate
    tab(ii,4)=mean(T(A(:)));      % mean dev. stress
    tab(ii,5)=std(T(A(:)));       % std dev. stress
    tab(ii,6)=mean(misfit(A(:))); % mean misfit
    tab(ii,7)=mean(taugii(A(:))); % mean grav. stress
    tab(ii,8)=mean(B(A(:)));      % mean bound. stress
    tab(ii,9)=circ_mean(2*shmax(A(:))*pi/180)/2*180/pi;
    tab(ii,10)=circ_confmean(2*shmax(A(:))*pi/180,0.05)/2*180/pi;      % std err of S-H-max azimuths
    tab(ii,11)=median(eta(A(:)));      % mean viscosity
end

file=fopen(strcat('summary-',model,'.csv'),'w');
s='region,points,strain,dev,devstd,misfit,grav,bound,shmax,shmaxste,visc';
fprintf(file,'%s\n',s);
t='%i,%i,%.2f,%.1f,%.1f,%.2f,%.1f,%.1f,%i,%.2f,%.1f\n';
fprintf(file,'%i,%i,%.2f,%.1f,%.1f,%.2f,%.1f,%.1f,%.0f,%.2f,%.1f\n',tab');
fclose(file);
% dlmwrite(strcat('summary-',model,'.csv'),tab,'-append','delimiter',',');      % write summary















