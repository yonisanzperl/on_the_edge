#!/bin/bash
#SBATCH --job-name=PertSub
#SBATCH --mail-type=END
#SBATCH --mail-user=
#SBATCH --mem-per-cpu=2G
#SBATCH --cpus-per-task=2
#SBATCH --array=1-176
#SBATCH --output=PertuSub_v1%A_%a.out
#SBATCH --error=PertuSub_v1%A_%a.err

#Load Matlab 2017a module
ml MATLAB

matlab -nojvm -nodisplay<<-EOF

s=str2num(getenv('SLURM_ARRAY_TASK_ID'))

%%  Read the empirical data 
load('SC_dtk68horn.mat')
%load('hcpunrelated100_REST_dkt68.mat')
load('DK68_empirical.mat')
load('hopf_freq_DKT68.mat')
C = Cnew;
C=C/max(max(C))*0.2;

N=68;
NPARCELLS = N;
NSUB=100;
NSUBSIM=100; 
Tmax=1200;
indexsub=1:NSUB;
Isubdiag = find(tril(ones(N),-1));

% 
% for nsub=indexsub
%     tsdata(:,:,nsub)=subject{1,nsub}.dkt68ts  ;
%     FCdata(nsub,:,:)=corrcoef(squeeze(tsdata(:,:,nsub))');
% end
% 
% FC_emp=squeeze(mean(FCdata,1));
% 
% FCemp2=FC_emp-FC_emp.*eye(N);
% GBCemp=mean(FCemp2,2);
% 
% Isubdiag = find(tril(ones(N),-1));
% 
% TR = 0.72
% %%%%%%%%%%%%%% filter for simulated signa "MEGstyle"
% 
% flp = 0.008;           % lowpass frequency of filter
% fhi = 0.08;           % highpass
% delt = TR;            % sampling interval
% k=2;                  % 2nd order butterworth filter
% fnq=1/(2*delt);       % Nyquist frequency
% Wn=[flp/fnq fhi/fnq]; % butterworth bandpass non-dimensional frequency
% [bfilt2,afilt2]=butter(k,Wn);   % construct the filter
% 
% 
% %%%%%%%%%%%%%%
% %%
% 
% kk=1;
% for nsub=1:NSUB
%     nsub
%     BOLDdata=(squeeze(tsdata(:,:,nsub)));
%     for seed=1:N
%         BOLDdata(seed,:)=BOLDdata(seed,:)-mean(BOLDdata(seed,:));
%         timeseriedata(seed,:) =filtfilt(bfilt2,afilt2,BOLDdata(seed,:));
%         Xanalytic = hilbert(demean(timeseriedata(seed,:)));
%         PhasesE(seed,:) = angle(Xanalytic);
%     end
%     FCdata_filt(nsub,:,:)=corrcoef(timeseriedata(:,:)');
%     
%     
% % Meta emp    
%     for i=1:NPARCELLS
%         gKOM=nansum(complex(cos(PhasesE),sin(PhasesE)))/NPARCELLS;
%         enstrophy1(i,:)=abs(gKOM);
%     end
%     
%     Metaemp(1,nsub)=nanstd(enstrophy1(:));
% 
%     
% 
%     ii2=1;
%     for t=1:15:Tmax-30
%         jj2=1;
%         cc=corrcoef((timeseriedata(:,t:t+30))');
%         for t2=1:15:Tmax-30
%             cc2=corrcoef((timeseriedata(:,t2:t2+30))');
%             ca=corrcoef(cc(Isubdiag),cc2(Isubdiag));
%             if jj2>ii2
%                 cotsamplingdata(kk)=ca(2);   %% this accumulate all elements of the FCD empirical
%                 kk=kk+1;
%             end
%             jj2=jj2+1;
%         end
%         ii2=ii2+1;
%     end
% end
% FC_empfilt = squeeze(nanmean(FCdata_filt,1)); 
% FCemp2filt=FC_empfilt-FC_empfilt.*eye(N);
% GBCempfilt=mean(FCemp2filt,2);
% Metaemp_av = nanmean(Metaemp);

%% the model
% fixed parameter
TR = 0.72;
% %%%%%%%%%%%%%% filter for simulated signa "MEGstyle"

flp = 0.008;           % lowpass frequency of filter
fhi = 0.08;           % highpass
delt = TR;            % sampling interval
k=2;                  % 2nd order butterworth filter
fnq=1/(2*delt);       % Nyquist frequency
Wn=[flp/fnq fhi/fnq]; % butterworth bandpass non-dimensional frequency
[bfilt2,afilt2]=butter(k,Wn);   % construct the filter


fcsimul=zeros(NSUB,NPARCELLS,NPARCELLS);



Tmax=1200;
omega1 = repmat(2*pi*f_diff',1,2); omega1(:,1) = -omega1(:,1);
dt=0.1*TR/2;
sig=0.0006;
dsig = sqrt(dt)*sig;


% Para el Hopf sin la diag
C = Cnew;
factor=max(max(C));
C=C/factor*0.2;


BETA=0:0.2:2;
%GG=0.:0.2:3;

GG=3.2:0.2:4;


[IG Ibeta]=ind2sub([length(GG),length(BETA)],s);
G=GG(IG);
beta=BETA(Ibeta);



wC = G*C;
sumC = repmat(sum(wC,2),1,2); % for sum Cij*xj
omega2=beta*ones(NPARCELLS,2);
omega2(:,1) = -omega2(:,1);
omega=omega1+omega2;
fcsimul=zeros(N,N,NSUBSIM);
kk=1;
for sub=1:NSUBSIM
    sub
    %% Hopf Simulation
    a=-0.02*ones(NPARCELLS,2);

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
        signal_filt(seed,:) =filtfilt(bfilt2,afilt2,ts(seed,:));
        Xanalytic = hilbert(demean(signal_filt(seed,:)));
        Phases(seed,:) = angle(Xanalytic);
    end
    
    fcsimul(:,:,sub)=corrcoef(signal_filt');
% Meta    
    for i=1:NPARCELLS
        gKOM=nansum(complex(cos(Phases(i,:)),sin(Phases(i,:))))/NPARCELLS;
        enstrophy1(i,:)=abs(gKOM);
    end
    
    Metasimul(1,sub)=nanstd(enstrophy1(:));
    
 % FCD   
    ii2=1;
    for t=1:15:Tmax-30
        jj2=1;
        cc=corrcoef((ts(:,t:t+30))');
        for t2=1:15:Tmax-30
            cc2=corrcoef((ts(:,t2:t2+30))');
            ca=corrcoef(cc(Isubdiag),cc2(Isubdiag));
            if jj2>ii2
                cotsamplingsim(kk)=ca(2);  %% FCD simulation
                kk=kk+1;
            end
            jj2=jj2+1;
        end
        ii2=ii2+1;
    end
    
end

%%
% Metric corr of FC
FC_simul=squeeze(mean(fcsimul,3));
cc=corrcoef(atanh(FC_empfilt(Isubdiag)),atanh(FC_simul(Isubdiag)));
FCfitt=cc(2); %% FC fitting
FCfitt2 = ssim(FC_simul, FC_empfilt);

FCsim2=FC_simul-FC_simul.*eye(N);
GBCsim=nanmean(FCsim2,2);
GBCfitt1=corr2(GBCempfilt,GBCsim);
GBCfitt2 = sqrt(nanmean((GBCempfilt-GBCsim).^2));

Metafitt = abs(nanmean(Metasimul)-Metaemp_av);

[hh pp FCDfitt]=kstest2(cotsamplingdata,cotsamplingsim);  %% FCD fitting

%%%
save(sprintf('PertSub68_op_%03d.mat',s),'GBCfitt1','GBCfitt2','FCfitt','FCfitt2', 'FCDfitt','Metafitt');

EOF