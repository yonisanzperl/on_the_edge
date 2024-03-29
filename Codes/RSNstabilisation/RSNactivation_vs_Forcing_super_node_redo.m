% #!/bin/bash
% #SBATCH --job-name=PertSub
% #SBATCH --mail-type=END
% #SBATCH --mail-user=
% #SBATCH --mem-per-cpu=2G
% #SBATCH --cpus-per-task=2
% #SBATCH --array=1-100
% #SBATCH --output=HopfInfoC_Sub%A_%a.out
% #SBATCH --error=HopfInfoC_Sub%A_%a.err
% 
% #Load Matlab 2017a module
% ml MATLAB
% 
% matlab -nojvm -nodisplay<<-EOF
% 
% s=str2num(getenv('SLURM_ARRAY_TASK_ID'))
% rng(s);

% the optimum G and Beta for Subcritics regimen
clear all
addpath('C:\Users\yonis\Documents\Laburo\PostDoc\PerturbationWPs\re_do_sub_super')
load('SC_dtk68horn.mat')
%load('hcpunrelated100_REST_dkt68.mat')
load('DK68_empirical_new.mat')
load('hopf_freq_DKT68.mat')

G= 0.4;
beta = 2.2; 

NPARCELLS=68;
N = NPARCELLS;
NSUBSIM=100;
NR=400;

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

%amp = 0:0.0005:0.01;

amp =0.02;

kick = repmat(ones(NPARCELLS,1),1,2);


wC = G*C;
sumC = repmat(sum(wC,2),1,2); % for sum Cij*xj
omega2=beta*ones(NPARCELLS,2);
omega2(:,1) = -omega2(:,1);
omega=omega1+omega2;
fcsimul=zeros(N,N,NSUBSIM);



for kk=1:NPARCELLS/2
    rng(kk);
    kick = zeros(NPARCELLS,1);
    kick(kk) = 1;
    kick(kk+34)=1;
    kick = repmat(kick,1,2);
    for sub=1:NSUBSIM
        sub
        
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
        tss=xs';
        for seed=1:NPARCELLS
            ts(seed,:)=detrend(ts(seed,:)-mean(ts(seed,:)));
            signal_filt(seed,:) =filtfilt(bfilt,afilt,ts(seed,:));
            Xanalytic = hilbert(demean(signal_filt(seed,:)));
            Phases(seed,:) = angle(Xanalytic);
        end
        KoP(sub)=mean(abs(sum(complex(cos(Phases),sin(Phases)),1))/N);
        FC_unp(:,:,sub)= corr(signal_filt');
        
        %%% Perturbation
        F0x = [amp,0];
        F0y = [0, amp];
        nn=0;
        
        
        for t=0:dt:((Tmax-1)*TR)
            suma = wC*z - sumC.*z; % sum(Cij*xi) - sum(Cij)*xj
            forcing = kick.*(F0x.*cos(omega_mean.*t)+F0y.*sin(omega_mean.*t));
            zz = z(:,end:-1:1); % flipped z, because (x.*x + y.*y)
            z = z + dt*(a.*z + zz.*omega - (z+omega2.*zz).*(z.*z+zz.*zz) + suma + forcing) + dsig*randn(NPARCELLS,2);
            if abs(mod(t,TR))<0.01
                nn=nn+1;
                xs(nn,:)=z(:,1)';
                forz(nn,:)=forcing(:,1);
                forz2(nn,:)=forcing(:,2);
            end
        end
        ts=xs';
        Rspatime1=zeros(NPARCELLS,Tmax);
        
        for seed=1:NPARCELLS
            ts(seed,:)=detrend(ts(seed,:)-mean(ts(seed,:)));
            signal_filt(seed,:) =filtfilt(bfilt,afilt,ts(seed,:));
            Xanalytic = hilbert(demean(signal_filt(seed,:)));
            Phases(seed,:) = angle(Xanalytic);
        end
        FC_pert(:,:,sub)= corr(signal_filt');
        KoPper(sub)=mean(abs(sum(complex(cos(Phases),sin(Phases)),1))/N);

    end
    infocapacity(kk)=std(KoPper-mean(KoP));
    susceptibility(kk)=mean(KoPper-mean(KoP));
    FC_pert_av(:,:,kk) =squeeze(nanmean(FC_pert,3));
    FC_unp_av(:,:,kk) =squeeze(nanmean(FC_unp,3));
    



save(sprintf('Wtrials_forz_sup_RSNact_%03d.mat',kk),'infocapacity','susceptibility','NSUBSIM','FC_pert_av','FC_unp_av');
end