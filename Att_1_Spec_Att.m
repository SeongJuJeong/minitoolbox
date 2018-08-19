clear all; close all;clc
addpath /Users/goodsjj/SMU/NorthTx/minitoolbox/
% cd /Users/goodsjj/SMU/NorthTx/Attenuation/FWB
Region0={'DFW'; 'Cleburne'; 'Venus'; 'Azle'; 'Irving'; 'NEIC'};
% Region0={'NEIC'};
delete *out NEIC/*out Azle/*out Irving/*out Venus/*out DFW/*out Cleburne/*out

for ii=1:length(Region0)
Region=Region0{ii};
for i=1:3
    KCMP =['SH';'SV';'UP'];
    CMP=KCMP(i,:);
    [t1 t2 t3 t4 t5 t6 t7 event t9 t10 t11]=....
    textread(['RawData/',Region,'_',num2str(CMP),'.lst'],...
    '%s %s %s %s %s %s %s %s %s %s %s',...
        'delimiter', '/');
% in MATLAB
clearvars -except i nn lst Region Region0 CMP t1 t2 t3 t4 t5 t6 t7 event t9 t10 t11 A B 
% in OCTAVE
%  clear -x i nn lst Region Region0 CMP

[a b]=size(event);
df0=0.01;   % reference frequency 

for n=1:a
    A{n} = readsac([t1{n},'/',t2{n},'/',t3{n},'/',t4{n},'/',...
        t5{n},'/',t6{n},'/',t7{n},'/',event{n},'/',t9{n},'/',t10{n},'/ins*',...
        num2str(CMP),'.acc.t.am']);
    B{n} = readsac([t1{n},'/',t2{n},'/',t3{n},'/',t4{n},'/',...
        t5{n},'/',t6{n},'/',t7{n},'/',event{n},'/',t9{n},'/',t10{n},'/ins*',...
        num2str(CMP),'.noise.dis.t.am']);
    
    %%%% signal spectrum
    spec = A{n}.DATA1;                  % signal acc spectrum (cm/s)
    df   = A{n}.DELTA;                  % (Hz) sampling in frequency
    distance = A{n}.DIST;               % (km) distance between station and event
    depth = A{n}.EVDP;                  % (km) depth
    Hypo0= sqrt(distance.^2+depth.^2);  % (km) hypocentral distance
    Eml  = A{n}.MAG;                    % local magnitude of event
    nyqf = A{n}.E;                      % (Hz) Nyquist frequency             
    sfreq = 0:df:nyqf;                  % (Hz) frequency

    %%%% noise spectrum
    nspec = B{n}.DATA1;                 % noise spectra (displacement)
    ndf   = B{n}.DELTA;                 % (Hz) sampling in noise spectra
    nfreq = 0:ndf:nyqf;                 % (Hz) noise frequency
    w=2.*pi.*nfreq;                     % omega=2*pi*f
    w2=w.^2;                            % omega^2
    nspec_acc=nspec.*w2';               % acceleration (cm/s)

    %%%% Smoothing
    specL=length(spec);                          % length of signal spectra
    nspecL=length(nspec_acc);                    % length of noise spectra
    smoothspec=smooth(spec,'moving');            % smoothing signal spectra with 5 point
    smoothnspec=smooth(nspec_acc,'moving');      % smoothing noise spectra with 5 point

    %%%% frequency coordinate using spline interpolation
    freq=0:df0:nyqf;                             % reference frequency with df=0.01 
    spec_spline0=spline(sfreq,smoothspec,freq);  % signal
    nspec_spline=spline(nfreq,smoothnspec,freq); % noise
    
    %%%% Calculate Signal-to-Noise Ratio
    s2n=20.*log10(spec_spline0./nspec_spline);
    s2n_sm=smooth(s2n,10,'moving');     % smoothing noise spectra with 10 point
    
    %%%% site correction
    if CMP == 'UP'
        spec_spline=spec_spline0;
    else
        if t7{n}(1:3) == 'NEI'
            spec_spline=spec_spline0;
        else
            [freq1 HV_avg_norm0]=textread(['/Users/goodsjj/SMU/NorthTx/Broad_band/RESULT/HV/',...
            t9{n},'/',t9{n},'_HV_site_ratio.dat'],'%f %f');
            HV_spline =spline(freq1,HV_avg_norm0,freq);  % signal
            spec_spline=spec_spline0./HV_spline;
        end
    end
    
    %%%% spacing of frequency band %%%%
    f0=logspace(log10(0.5),log10(50),40);
    if nyqf > 99
        f=f0;
    elseif nyqf == 50
        f=f0(1:34);
    elseif nyqf == 20
        f=f0(1:26);
    end
    
    %%%% Choose amplitude at 40 frequencies
    iff=round(f./df0) +1;
    Att=spec_spline(iff);
    s2n_chosen=real(s2n_sm(iff));
   
    %%%% save result
    Hypo=Hypo0.*ones(length(f),1);                          % Hypo dist. (km)
%     V0=A{n}.USER9;
    [Vp]=GeoAzle2(depth,Region); % P-velocity (km/s)
    if CMP == 'UP';   V0=Vp;  else;    V0=Vp/1.8;  end
    V=V0.*ones(length(f),1);                                % velocity (km/s)
    for i=1:length(f);  llst(i,:)=event{n};  end            % event date
    f=reshape(f,length(f),1);                               % chosen frequency
    Att=reshape(Att,length(Att),1);                         % chosen spectral amp.
    s2n_chosen=reshape(s2n_chosen,length(s2n_chosen),1);    % chosen snr
    result=[llst,char(' '*char(ones(length(f),1))),...
        num2str(V,'%.1f'),char(' '*char(ones(length(f),1))),num2str(Hypo,'%.1f'),...
        char(' '*char(ones(length(f),1))),num2str(f,'%.2f'),...
        char(' '*char(ones(length(f),1))),num2str(Att,'%.1e'),...
        char(' '*char(ones(length(f),1))),num2str(s2n_chosen,'%.2f')];
    dlmwrite([Region,'/',event{n},'_',num2str(CMP),'_Att.out'],...
        result,'delimiter','','-append')
    % In OCTAVE
%     dlmwrite([num2str(lst{nn}),'_',num2str(CMP),'_mw.out'],...
%     result,"delimiter",r,"-append")
clear llst
end
end
end

% !cat *SH_Att.out > SH_f.out
% !cat *SV_Att.out > SV_f.out

!cat */*SH_Att.out */*SV_Att.out > SS_f.out
!cat */*UP_Att.out > UP_f.out

Att_2_Castro_Oth_att
