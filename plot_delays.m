clc; close all; clear all;

obsfile = 'dt_tst_utm.nc';
obsfile_geo = '../MATLAB/matgsdf-master/dt_tst.nc';
modfile = 'delays.out';
dummy = -999.9;

obs = ncread(obsfile,'dtp');
stla = ncread(obsfile_geo,'MeasPosX');
stlo = ncread(obsfile_geo,'MeasPosY');
evlo = ncread(obsfile_geo,'EventPosY');
evla = ncread(obsfile_geo,'EventPosX');
staco = ncread(obsfile,'StatComb');
estaco = ncread(obsfile,'EventStatComb');
periods =  ncread(obsfile,'Periods');
dims = size(obs);
nrays = dims(1);
nevents = dims(2);
nperiods = length(periods);

for p=1:nperiods
    
    if p<4
        mod = dlmread(modfile,'',[1+(p-1)*nrays*nevents, 0, p*nrays*nevents, 4]);
        mod(:,1:3) = mod(:,1:3) + 1;
        
        [statdist,~]=distance(stla(mod(mod(:,1)~=dummy & mod(:,5)~=dummy,2)), stlo(mod(mod(:,1)~=dummy & mod(:,5)~=dummy,2)), stla(mod(mod(:,1)~=dummy & mod(:,5)~=dummy,3)), stlo(mod(mod(:,1)~=dummy & mod(:,5)~=dummy,3)),'degrees');
        statdist = 6371 .* deg2rad(statdist);
        
        [ddist1, ~] = distance(evla(mod(mod(:,1)~=dummy & mod(:,5)~=dummy,1)), evlo(mod(mod(:,1)~=dummy & mod(:,5)~=dummy,1)), stla(mod(mod(:,1)~=dummy & mod(:,5)~=dummy,2)), stlo(mod(mod(:,1)~=dummy & mod(:,5)~=dummy,2)),'degrees');
        [ddist2, ~] = distance(evla(mod(mod(:,1)~=dummy & mod(:,5)~=dummy,1)), evlo(mod(mod(:,1)~=dummy & mod(:,5)~=dummy,1)), stla(mod(mod(:,1)~=dummy & mod(:,5)~=dummy,3)), stlo(mod(mod(:,1)~=dummy & mod(:,5)~=dummy,3)),'degrees');
        ddist1 = 6371 .* deg2rad(ddist1);
        ddist2 = 6371 .* deg2rad(ddist2);
        ddist2 = ddist1 - ddist2;
        
        mod = reshape(mod(:,end),nevents,nrays)';
        mod = mod(:);
        mod(mod==dummy) = [];
        
        obs_p = obs(:,:,p);
        obs_lst = obs_p';
        obs_lst = obs_lst(:);
        obs_lst(isnan(obs_lst)) = [];
        obs_p = obs_p(:);
        obs_p(isnan(obs_p)) = [];
        
        figure; set(gcf,'position',[1,1,960,918]);
        subplot(3,1,1); hold on;
        plot(obs_p, mod, '+k')
        plot([-100 100], [-100 100], 'r', 'LineWidth', 2)
        title(['Period: ' num2str(periods(p)) ' s']);
        xlim([-200 200])
        ylim([-100 100])
        xlabel('Observed \delta\tau_{ph} [sec]')
        ylabel('Modelled \delta\tau_{ph} [sec]')
        
    else
        
        [statdist,~]=distance(stla(mod(mod(:,1)~=dummy & mod(:,5)~=dummy,2)), stlo(mod(mod(:,1)~=dummy & mod(:,5)~=dummy,2)), stla(mod(mod(:,1)~=dummy & mod(:,5)~=dummy,3)), stlo(mod(mod(:,1)~=dummy & mod(:,5)~=dummy,3)),'degrees');
        statdist = 6371 .* deg2rad(statdist);
        
        [ddist1, ~] = distance(evla(mod(mod(:,1)~=dummy & mod(:,5)~=dummy,1)), evlo(mod(mod(:,1)~=dummy & mod(:,5)~=dummy,1)), stla(mod(mod(:,1)~=dummy & mod(:,5)~=dummy,2)), stlo(mod(mod(:,1)~=dummy & mod(:,5)~=dummy,2)),'degrees');
        [ddist2, ~] = distance(evla(mod(mod(:,1)~=dummy & mod(:,5)~=dummy,1)), evlo(mod(mod(:,1)~=dummy & mod(:,5)~=dummy,1)), stla(mod(mod(:,1)~=dummy & mod(:,5)~=dummy,3)), stlo(mod(mod(:,1)~=dummy & mod(:,5)~=dummy,3)),'degrees');
        ddist1 = 6371 .* deg2rad(ddist1);
        ddist2 = 6371 .* deg2rad(ddist2);
        ddist2 = ddist1 - ddist2;
        
        obs_p = obs(:,:,p);
        obs_lst = obs_p';
        obs_lst = obs_lst(:);
        obs_lst(isnan(obs_lst)) = [];
        obs_p = obs_p(:);
        obs_p(isnan(obs_p)) = [];
        
        figure; set(gcf,'position',[1,1,960,918]);
        subplot(3,1,1);
        plot(0, 0, '+k')
        
    end
    
    subplot(3,1,2)
    plot(statdist, obs_lst, '+k')
    xlabel('Great circle distance [km]')
    ylabel('Observed \delta\tau_{ph} [sec]')
    ylim([-200 200])
    
    subplot(3,1,3)
    plot(ddist2, obs_lst, '+k')
    xlabel('\Delta dist_{epi} [km]')
    ylabel('Observed \delta\tau_{ph} [sec]')
    ylim([-200 200])
    
    
    
    print(num2str(periods(p)),'-dpng')
    close all
    
end