% Source parameter calculation in time domain
% Written by Jeong,S.J. 
% reference from Urbancic et al., BSSA 1996
%=========================================================================== 

function [E0 ro ohm01 fc2 scaledsig s_moment1 amw1 Tsnoke] = time_moment...
    (V,r,fs,dist,nyqf,K,SD,SV,f1,f2,fc1,CMP,freq)        
% Constant 
rho=2.68;                                   % g/cm^3               
beta=V*10.0^5;                              % km/sec -> cm/sec     
dt=1/fs;                                    % sampling period
fn=nyqf;                                    % nyquist frequency
distance=dist*10.0^5;                   % km -> cm 
[gs]=geo_att(distance,CMP,150);     	    % geometrical spreading
    
% kappa: kq
kq_P=0.00018055;
kq_S=0.00034348;

% Quality factor from kq
Qp=1/(kq_P*V);
Qs=1/(kq_S*V);

% local kappa: kq+ks
kappa_local=0.027637;

% Q value 
if CMP == 'UP'
    q=exp(-pi*distance*fc1/(beta*Qp));
    q0=Qp;
else
    if dist < 25
        q=exp(-pi*kappa_local*fc1);
        q0=distance/(kappa_local*beta);
    else
        q=(exp(-pi*distance*fc1/(beta*Qs)));
        q0=Qs;
    end
end

site_term=(exp(-pi*fc1*kappa));

% filter
[b,a]=butter(3,[f1/fn f2/fn]);
% SD1=filtfilt(b,a,SD)*gs/q;
% SV1=filtfilt(b,a,SV)*gs/q;
SD1=filtfilt(b,a,SD)*gs;
SV1=filtfilt(b,a,SV)*gs;

% energy flux
SD2=0.0;  SV2=0.0;
for i=1:length(SD1)
    SD2=SD2+SD1(i)^2*dt;
end
for i=1:length(SV1)
    SV2=SV2+SV1(i)^2*dt;
end
SD2=2.0*SD2;					% cm^2*s
SV2=2.0*SV2;					% cm^2/s, energy flux

% source parameters
ohm01=sqrt(4*(SD2^(3/2)*SV2^(-1/2)));		
fc2=1/(2*pi)*sqrt(SV2/SD2);                 % Hz, corner frequency
s_moment1=(4*pi*rho*beta^3*10.0^5*ohm01)/r;	% dyne-cm, seismic moment
amw1=2.0/3.0*log10(s_moment1)-10.7; 		% moment magnitude
E0=4*pi*rho*beta*distance^2*SV2;            % g*cm^2/s^2, seismic energy
ro=(K*beta)/(2*pi*fc2);                     % cm, source radius
dsig1=(7*s_moment1)/(16*ro^3);              % dyne/cm^2, stress drop  
scaledsig=dsig1/1.0e6;                      % bars, stress drop
if scaledsig > 1000
    scaledsig=0;
end
                                                               
% Tsnoke=ohm01./(1.0+(freq./fc2).^2)...
%    .*exp(-pi.*distance.*freq./(beta.*q0));
Tsnoke=ohm01./(1.0+(freq./fc2).^2);
