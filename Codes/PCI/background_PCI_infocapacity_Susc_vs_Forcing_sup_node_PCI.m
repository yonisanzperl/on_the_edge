% the optimum G and Beta for Subcritics regimen



load('SC_dtk68horn.mat')
%load('hcpunrelated100_REST_dkt68.mat')
load('DK68_empirical.mat')
load('hopf_freq_DKT68_bad.mat')
addpath('C:\Users\yonis\Documents\Laburo\PostDoc\utils')

G= 0.4;
beta = 2.2;

NPARCELLS=68;
N = NPARCELLS;
NSUBSIM=30;

C = Cnew;
Isubdiag = find(tril(ones(NPARCELLS),-1));
%%%
% Parameters of the data
TR=0.72;  % Repetition Time (seconds)

% Bandpass filter settings
fnq=1/(2*TR);                 % Nyquist frequency
flp = 0.008;                    % lowpass frequency of filter (Hz)
fhi = 0.08;                    % highpass
Wn=[flp/fnq fhi/fnq];         % butterworth bandpass non-dimensional frequency
k=2;                          % 2nd order butterworth filter
[bfilt,afilt]=butter(k,Wn);   % construct the filter
Isubdiag = find(tril(ones(NPARCELLS),-1));

% Parameters HOPF
Tmax=1200;
omega1 = repmat(2*pi*f_diff',1,2); omega1(:,1) = -omega1(:,1);
omega_mean = mean(omega1(:,2),1);
dt=0.1*TR/2;
sig=0.01;
dsig = sqrt(dt)*sig;

%%

lam_mean_spatime_enstrophy=zeros(NPARCELLS,Tmax);
ensspasub=zeros(NSUBSIM,NPARCELLS);
ensspasub1=zeros(NSUBSIM,NPARCELLS);




factor=max(max(C));
C=C/factor*0.2;
% defines the forcing



wC = G*C;
sumC = repmat(sum(wC,2),1,2); % for sum Cij*xj
omega2=beta*ones(NPARCELLS,2);
omega2(:,1) = -omega2(:,1);
omega=omega1+omega2;
fcsimul=zeros(N,N,NSUBSIM);

for sub=1:NSUBSIM
    
    sub
    clear ts xs signal_filt Xanalytic Phases
    %% Hopf Simulation
    a=1.3*ones(NPARCELLS,2);
    xs=zeros(Tmax,NPARCELLS);
    %number of iterations, 100 willk�hrlich, weil reicht in diesem Fall
    z = 0.1*ones(NPARCELLS,2); % --> x = z(:,1), y = z(:,2)
    nn=0;
    % discard first 2000 time steps
    for t=0:dt:2000
        suma = wC*z - sumC.*z; % sum(Cij*xi) - sum(Cij)*xj
        zz = z(:,end:-1:1); % flipped z, because (x.*x + y.*y)
        z = z + dt*(a.*z + zz.*omega - (z+omega2.*zz).*(z.*z+zz.*zz) + suma) + dsig*randn(NPARCELLS,2);
    end
    % actual modeling (x=BOLD signal (Interpretation), y some other oscillation)
    for t=0:dt:((Tmax-1)*TR)
        suma = wC*z - sumC.*z; % sum(Cij*xi) - sum(Cij)*xj
        zz = z(:,end:-1:1); % flipped z, because (x.*x + y.*y)
        z = z + dt*(a.*z + zz.*omega - (z+omega2.*zz).*(z.*z+zz.*zz) + suma) + dsig*randn(NPARCELLS,2);
        if abs(mod(t,TR))<0.01
            nn=nn+1;
            xs(nn,:)=z(:,1)';
        end
    end
    ts=xs';
    
    for seed=1:NPARCELLS
        ts(seed,:)=detrend(ts(seed,:)-mean(ts(seed,:)));
        signal_filt(seed,:) =filtfilt(bfilt,afilt,ts(seed,:));
        Xanalytic = hilbert(demean(signal_filt(seed,:)));
        Phases(seed,:) = angle(Xanalytic);
    end
    

    
    
    % PCI exactly after perturbation
    ts_z = zscore(ts(:,1:200)); % zscore to all ts together and short
    ts_z(find(ts_z<=2))=0;
    ts_z(find(ts_z>2))= 1;
    [CC, ~, ~] = calc_lz_complexity(ts_z(:), 'exhaustive', true);
    PCI_sub_short(sub,1) = CC
    
    
    %susceptibility only after perturbation (Phases only includes
    %afteer pertubation signal)
    
end


PCI_short= nanmean(PCI_sub_short,1);

%save(sprintf('Wtrials_forz_sub_nodePCIredo_%03d.mat',s),'infocapacity','susceptibility','PCI_short');

