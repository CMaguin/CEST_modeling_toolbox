function rfmt=RF_MT(T2c,w1,dw,lineshape,cutoff)

if nargin<5; cutoff=gamma_2pi*11.7; end % by default the Superlorentzian cutoff is +-1ppm for 11.7T

if strcmp(lineshape,'SuperLorentzian') %SuperLorentzian
    
    for i=1:numel(w1)
        rfmt(:,i)=superLorentzian_extrap_DG(T2c,w1(i),dw, cutoff);
    end

    rfmt=rfmt';

elseif strcmp(lineshape,'Gaussian') %Gaussian
    rfmt=w1.^2*T2c.*sqrt(pi/2).*exp(-(dw.*T2c).^2./2); 
elseif strcmp(lineshape,'Lorentzian') %Lorentzian
    rfmt=w1.^2*T2c./(1+(dw.*T2c).^2);
elseif strcmp(lineshape,'SL_CorrectedDipolarOrder') %Adding contribution made by dipolar order
    for i=1:numel(w1)
        rfmt(:,i)=superLorentzian_extrap_DG(T2c,w1(i),dw, cutoff);
    end

    rfmt=rfmt';
    
    T1D=6.1*10^(-3);%5,6ms in GM, 6.1 ms in WM
    D=sqrt(1/15)/T2c; %from morrison et al 1995
    rfmt=rfmt ./(1+T1D.*rfmt.*(dw/D).^2) ;%steady state solution
else
    error('Unknown MT-lineshape - choose SuperLorentzian, Lorentzian or Gaussian');
end;


    