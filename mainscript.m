% A simplified "Fast" Monte Carlo generator for understanding
% fiducial phase space, efficiency, kinematic effects at threshold etc.
%
% Example use cases:
%
%    - Test difference between e.g.
%      |y_system| < 0.9 VS |eta_pi| < 0.9 & pT > 0.2 GeV phase spaces
%
%    - Test how the invariant mass distribution changes, when the
%      pseudorapidity cut of final states is changed
%
%    - Test difference between Lorentz frames with single limits Pt -> 0,
%      Y -> 0 or double limit Pt and Y -> 0 (all frames should coincide)
% 
% If you port to e.g. C++, use quality C++11 random numbers / ROOT TRandom3
%
% 4-momentum convention is p = [px,py,pz,E] = [p(1),p(2),p(3),p(4)]
%
% mikael.mieskolainen@cern.ch, 13/07/2018
clear; close all;

rng('default');   % Random numbers

mpi = 0.139570;   % Charged pion mass
mK  = 0.493677;   % Charged kaon mass

% Minkowski metric
g = [-1  0  0  0;
     0  -1  0  0;
     0   0 -1  0
     0   0  0  1];
 
% Load kinematic functions
kinfunctions;


%% Generator parameters

PTMODE   = 2;     % 1 for flat in pt^2, 2 for exponential in pt^2, 3 for flat in pt
MASSMODE = 3;     % 1 for flat in m^2,  2 for exponential in m^2,  3 for flat in m

lambda = 1 / 0.2; % System pt parameter (if PTMODE = 2)
kappa  = 1 / 1.2; % System mass parameter (if MASSMODE = 2)

mdec = [mpi mpi]; % Decay daughter masses (2 x 1)


%% Event generation sampling boundaries (must be set)

% System pt (GeV)
limits.ptmin = 0.0;
limits.ptmax = 3.0;

% System rapidity
limits.ymin = -0.9;
limits.ymax =  0.9;

% System invariant mass (GeV)
limits.mmin = sum(mdec); % sum of decay daughter masses
limits.mmax = 2.5;


%% Fiducial acceptance cuts

fcuts_on = true;    % Cuts on/off (PLAY WITH THIS FIRST!)
etamax   = 0.9;     % Absolute pseudorapidity of final state particles
ptmin    = 0.175;   % Minimum pt of final state particles


%% Event loop

events = 1e4; % Number of events

% Observables we save in these vectors
ptvals = zeros(events,1);
mvals  = zeros(events,1);
yvals  = zeros(events,1);

piptvals  = zeros(events,1);
pietavals = zeros(events,1);
piphivals = zeros(events,1);

eta12vals      = zeros(events,2);
costhphivals   = zeros(events,2);
costhphivalsGJ = zeros(events,2);
costhphivalsHE = zeros(events,2);
costhphivalsCS = zeros(events,2);
costhphivalsBE = zeros(events,2);

generated = 0;

k = 1;
while (k < events + 1)
    
    % Generate event
    [p,p1,p2] = generator(PTMODE, MASSMODE, lambda, kappa, limits, mdec);
    
    % Apply fiducial cuts
    all_fid_ok = (f_pt(p1) > ptmin) && (f_pt(p2) > ptmin) && ...
                 (abs(f_eta(p1)) < etamax) && (abs(f_eta(p2)) < etamax);
    
    if (all_fid_ok || (fcuts_on == false))
        
        % System observables
        ptvals(k) = f_pt(p);
        mvals(k)  = f_m(p);
        yvals(k)  = f_rap(p);
        
        % Decay daughter 1D observables
        piptvals(k)  = f_pt(p1);
        pietavals(k) = f_eta(p1);
        piphivals(k) = f_phi(p1);
        
        % Decay daughter 2D observables
        eta12vals(k,:) = [f_eta(p1) f_eta(p2)];
        
        % Rest frame (no rotations here)
        sign = -1;
        p1rf = boostroutine(p, p1, sign);
        p2rf = boostroutine(p, p2, sign);
        costhphivals(k,:)   = [cos(f_theta(p1rf)) f_phi(p1rf)];
        
        % GJ-frame
        pf = {p1, p2}; direction = 1; sqrts = 7000;
        pfout = GJframe(pf, direction, sqrts);
        costhphivalsGJ(k,:) = [cos(f_theta(pfout{1})) f_phi(pfout{1})];
        
        % HE-frame
        pfout = HEframe(pf, direction, sqrts);
        costhphivalsHE(k,:) = [cos(f_theta(pfout{1})) f_phi(pfout{1})];
        
        % CS-frame
        pfout = CSframe(pf, direction, sqrts);
        costhphivalsCS(k,:) = [cos(f_theta(pfout{1})) f_phi(pfout{1})];
        
        % BE-frame
        pfout = BEframe(pf, direction, sqrts);
        costhphivalsBE(k,:) = [cos(f_theta(pfout{1})) f_phi(pfout{1})];
        
        k = k + 1;
        if (mod(k, 1000) == 0)
            fprintf('Event %d/%d \n', k, events);
        end
    end
    generated = generated + 1;
end

fprintf('Fiducial cuts: etamax = %0.2f, ptmin = %0.2f GeV \n', etamax, ptmin);
fprintf('Fiducial Accepted / Total Generated = %0.3f \n', (k-1) / generated);
fprintf('Mean system pt = %0.2f GeV \n', mean(ptvals));


%% Plot histograms
% [replace hist() below with proper histogram functions from
%  https://github.com/mieskolainen/matlabcodes . hist() ok for quick tests.]
close all;

figure;
subplot(3,3,1);
hist(ptvals, 100); axis square; xlabel('system $P_t$ (GeV)','interpreter','latex'); axis([0 2 0 inf]);

subplot(3,3,2);
hist(mvals, 100); axis square; xlabel('system $M$ (GeV)','interpreter','latex'); axis([0 limits.mmax*1.1 0 inf]);

subplot(3,3,3);
hist(yvals, 100); axis square; xlabel('system rapidity $y$','interpreter','latex'); axis([limits.ymin*1.5 limits.ymax*1.5 0 inf]);

subplot(3,3,4);
hist(piptvals, 100); axis square; xlabel('daughter $p_t$ (GeV)','interpreter','latex'); axis([0 2.5 0 inf]);

subplot(3,3,5);
hist(pietavals, 100); axis square; xlabel('daughter $\eta$','interpreter','latex'); axis([-3.0 3.0 0 inf]);

subplot(3,3,6);
hist(piphivals, 100); axis square; xlabel('daughter $\phi$ (rad)','interpreter','latex'); axis([-pi-0.5 pi+0.5 0 inf]);

subplot(3,3,7);
[X,bins] = hist3(eta12vals, [50 50]); imagesc(bins{1}, bins{2}, X'); axis square;
xlabel('daughter$^+$ $\eta$','interpreter','latex');
ylabel('daughter$^-$ $\eta$','interpreter','latex');
set(gca,'yDir','normal'); axis([-2.0 2.0 -2.0 2.0]);
colormap(hot);

% Print out pdf
outputstr = sprintf('./pdf/output_kinematics_fiducial_%d.pdf', fcuts_on);
eval(sprintf('print -dpdf %s', outputstr));
system(sprintf('pdfcrop --margins 10 %s %s', outputstr, outputstr));

figure;

subplot(1,5,1);
[X,bins] = hist3(costhphivals, [50 50]); imagesc(bins{1}, bins{2}, X'); axis square;
xlabel('$\cos \theta$','interpreter','latex','fontsize',8);
ylabel('$\phi$ (rad)','interpreter','latex');
title(sprintf('RF, $S = %0.2f$', shannonentropy(X)),'interpreter','latex','fontsize',8);
set(gca,'yDir','normal'); axis([-1.0 1.0 -pi pi]);

subplot(1,5,2);
[XGJ,bins] = hist3(costhphivalsGJ, [50 50]); imagesc(bins{1}, bins{2}, XGJ'); axis square;
xlabel('$\cos \theta$','interpreter','latex','fontsize',8);
ylabel('$\phi$ (rad)','interpreter','latex');
title(sprintf('GJP, $S = %0.2f$', shannonentropy(XGJ)),'interpreter','latex','fontsize',8);
set(gca,'yDir','normal'); axis([-1.0 1.0 -pi pi]); colormap('hot');

subplot(1,5,3);
[XHE,bins] = hist3(costhphivalsHE, [50 50]); imagesc(bins{1}, bins{2}, XHE'); axis square;
xlabel('$\cos \theta$','interpreter','latex','fontsize',8);
ylabel('$\phi$ (rad)','interpreter','latex');
title(sprintf('Helicity, $S = %0.2f$', shannonentropy(XHE)),'interpreter','latex','fontsize',8);
set(gca,'yDir','normal'); axis([-1.0 1.0 -pi pi]); colormap('hot');

subplot(1,5,4);
[XCS,bins] = hist3(costhphivalsCS, [50 50]); imagesc(bins{1}, bins{2}, XCS'); axis square;
xlabel('$\cos \theta$','interpreter','latex','fontsize',8);
ylabel('$\phi$ (rad)','interpreter','latex');
title(sprintf('Collins-Soper, $S = %0.2f$', shannonentropy(XCS)),'interpreter','latex','fontsize',8);
set(gca,'yDir','normal'); axis([-1.0 1.0 -pi pi]); colormap('hot');

subplot(1,5,5);
[XBE,bins] = hist3(costhphivalsBE, [50 50]); imagesc(bins{1}, bins{2}, XBE'); axis square;
xlabel('$\cos \theta$','interpreter','latex','fontsize',8);
ylabel('$\phi$ (rad)','interpreter','latex');
title(sprintf('Anti-Helicity, $S = %0.2f$', shannonentropy(XBE)),'interpreter','latex','fontsize',8);
set(gca,'yDir','normal'); axis([-1.0 1.0 -pi pi]); colormap('hot');

% Print out pdf
outputstr = sprintf('./pdf/output_frames_fiducial_%d.pdf', fcuts_on);
eval(sprintf('print -dpdf %s', outputstr));
system(sprintf('pdfcrop --margins 10 %s %s', outputstr, outputstr));

