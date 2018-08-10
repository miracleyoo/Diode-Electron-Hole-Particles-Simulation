function [scatGaAs_hole,Gmh,Gml]=make_GaAs_hole_scatTable_v2(T,de,Vmax,cdop)

%--------Electron Energy Steps for Formulas/Graphs-------------------------
delt_Ek=de;                     %Electron Energy Step (default 0.002);
Ek_pts=Vmax/de;                 %Number of Sample Points for Electron Energy   
eV_axis=(1:Ek_pts)*delt_Ek;
scat_h(Ek_pts,12)=0;              %Matrix holding scattering rates for Gamma Valley
scat_l(Ek_pts,12)=0;              %Matrix holding scattering rates for L Valley

%---------General Constants------------------------------------------------
q=1.60219e-19;                  %Charge of Electron
kb=1.38066e-23;                 %Boltzmann's Constant
h=1.05459e-34;                  %Planck's Constant (/2pi)
hwo=0.03536;                    %In eV -> Need 0.03536*q for calculations
hwoq=hwo*q;
wo=hwoq/h;  
emR=9.10953e-31;                %Mass of Electron at Rest
eps_o=8.85419e-12;              %Vacuum Permittivity
eps_s=12.9;                     %for GaAs,es
eps_inf=10.92;                  %Optical Permittivity for GaAs
Ni=cdop;                        %m^-3, ionized density
m=emR;%free-electron mass
ml=0.082*m;%light hole
mh=0.45*m;%heavy hole
hw0=0.035;%eV
w0=hw0*1.6e-19/h;
adp=7*q;                        %Acoustic Deformation Potential
rho=5360;%kg/m^3
s=3860;%m/s
%inverse band mass parameters
A=-6.98;
B=-4.5;
C=6.2;

%-------acoustic phonon scattering-----------------------------------------
for i=1:Ek_pts
    ei=delt_Ek*i;
    ef=ei;%consider as elastic
    gamma=adp*adp*kb*T/((2*pi)^2*h*rho*s*s);
    a=h*h*abs(A)/(2*emR);
    %k_theta=;
    g=sqrt((B/A)^2+(C/A)^2/3);%choose the max g in heavy hole
    %P=gamma*ei^0.5*0.25*(1+3*cos(k_theta)^2)/(a*(1-g))^1.5;
    P=gamma*ei^0.5*0.5*1/(a*(1-g))^1.5;%use the max overlap factor
    scat_h(i,11)=P/1e9;%h-h and 
    scat_l(i,12)=P/1e9;%l-h
    g=B/A;
    P=gamma*ei^0.5*0.5*1/(a*(1+g))^1.5;%use the max overlap factor
    scat_l(i,11)=P/1e9;%l-l and 
    scat_h(i,12)=P/1e9;%h-l
end

%nonoplar optical phonon scattering
DK2=1.58e22; %ev^2/m^2
for i=1:Ek_pts
    %final is heavy,abs
    ei=delt_Ek*i;
    ef=ei+hw0;
    N0 = 1 / (exp( (q*hw0)/(kb*T) ) - 1);
    B0=h^2*DK2/(2*rho*hw0);
    gamma=2*pi/h*B0*N0;
    a=h*h*abs(A)/(2*emR);
    g=sqrt((B/A)^2+(C/A)^2/3);%choose the max g in heavy hole
    P=gamma/(2*pi)^3*ef^0.5*0.5/(a*(1-g))^1.5;
    scat_h(i,7)=P/1e27;%h-h,absorbsion
    scat_l(i,7)=P/1e27;%l-h,absorbsion
    %final is light, absorbsion
    ei=delt_Ek*i;
    ef=ei+hw0;
    N0 = 1 / (exp( (q*hw0)/(kb*T) ) - 1);
    B0=h^2*DK2/(2*rho*hw0);
    gamma=2*pi/h*B0*N0;
    a=h*h*abs(A)/(2*emR);
    g=B/A;
    P=gamma/(2*pi)^3*ef^0.5*0.5/(a*(1+g))^1.5;
    scat_h(i,9)=P/1e27;%h-l,absorbsion
    scat_l(i,9)=P/1e27;%l-l,absorbsion
    %final is heavy, emission
    ei=delt_Ek*i;
    if (ei-hw0)<=de
        P=0;
        scat_h(i,8)=P;
        scat_l(i,8)=P;
    else
        ef=ei-hw0;
        N0 = 1 / (exp( (q*hw0)/(kb*T) ) - 1);
        B0=h^2*DK2/(2*rho*hw0);
        gamma=2*pi/h*B0*(N0+1);
        a=h*h*abs(A)/(2*emR);
        g=sqrt((B/A)^2+(C/A)^2/3);%choose the max g in heavy hole
        P=gamma/(2*pi)^3*ef^0.5*0.5/(a*(1-g))^1.5;
        scat_h(i,8)=P/1e27;%h-h,emisson
        scat_l(i,8)=P/1e27;%l-h,emisson
    end
    %final is light, emission
    ei=delt_Ek*i;
    if (ei-hw0)<=de
        P=0;
        scat_h(i,10)=P;
        scat_l(i,10)=P;
    else
        ef=ei-hw0;
        N0 = 1 / (exp( (q*hw0)/(kb*T) ) - 1);
        B0=h^2*DK2/(2*rho*hw0);
        gamma=2*pi/h*B0*(N0+1);
        a=h*h*abs(A)/(2*emR);
        g=B/A;
        P=gamma/(2*pi)^3*ef^0.5*0.5/(a*(1+g))^1.5;
        scat_h(i,10)=P/1e27;%h-l,emisson
        scat_l(i,10)=P/1e27;%l-l,emisson
    end
end

%--------Impurity Scattering-----------------------------------------------
for i=1:Ek_pts
    %heavy holes to light holes---------------------------------
    ki=sqrt(2*mh*delt_Ek*i*q/(h*h));%related to E---------------------
    kf=sqrt(ml/mh)*ki;
    beta=sqrt(Ni*q*q/(kb*T*eps_o*eps_s));
    F=(beta*beta+ki*ki+kf*kf)/(ki*kf)*log((beta*beta+(ki+kf)^2)/(beta*beta+(ki-kf)^2))-4;
    P=3*q^4*Ni*ml*F/(32*pi*h^3*eps_o^2*eps_s^2*ki^2*kf);
    scat_h(i,1)=P;
    %light to heavy holes----------------------------------
    ki=sqrt(2*ml*delt_Ek*i*q/(h*h));
    kf=sqrt(mh/ml)*ki;
    F=(beta*beta+ki*ki+kf*kf)/(ki*kf)*log((beta*beta+(ki+kf)^2)/(beta*beta+(ki-kf)^2))-4;
    P=3*q^4*Ni*mh*F/(32*pi*h^3*eps_o^2*eps_s^2*ki^2*kf);
    scat_l(i,1)=P;
    %heavy to heavy-----------------------
    ki=sqrt(2*mh*delt_Ek*i*q/(h*h));
    kf=ki;
    F=(beta*beta+2*ki*ki)/(ki*ki)*log((beta*beta)/(beta*beta+4*ki*ki))+4/3*(3*beta^4+12*beta^2*ki^2+8*ki^4)/(beta^2*(beta^2+4*ki^2));
    P=3*q^4*Ni*mh*F/(32*pi*h^3*eps_o^2*eps_s^2*ki^2*kf);
    scat_h(i,2)=P;
    %light to light--------------------
    ki=sqrt(2*ml*delt_Ek*i*q/(h*h));
    kf=ki;
    F=(beta*beta+2*ki*ki)/(ki*ki)*log((beta*beta)/(beta*beta+4*ki*ki))+4/3*(3*beta^4+12*beta^2*ki^2+8*ki^4)/(beta^2*(beta^2+4*ki^2));
    P=3*q^4*Ni*ml*F/(32*pi*h^3*eps_o^2*eps_s^2*ki^2*kf);
    scat_l(i,2)=P;
end

%optical phonon scattering
for i=1:Ek_pts
%heavy holes to heavy holes,aborbtion---------------------------------
    ep = 1 / ( (1 / (eps_inf*eps_o)) - (1 / (eps_s*eps_o)) );
    N0 = 1 / (exp( (q*hw0)/(kb*T) ) - 1);
    c_const =(q/h)^3 * (hw0 * mh) / (4 * pi * ep);
    ab_const = c_const * N0;
    ei=delt_Ek*i;
    ef=ei+hw0;
    ki=sqrt(2*mh*delt_Ek*i*q/(h*h));
    kf=sqrt(2*mh*ef*q/(h*h));
    %kf=sqrt(ml/mh)*(ki^2+2*mh*w0/h);
    fai = log((ki + kf) / (max(ki, kf) - min(ki, kf)));        
    theta = (ki * ki + kf * kf) / (2 * ki * kf); 
    hi = (1 + 3 * theta * (theta - 1/fai)) / 4; 
    P=ab_const * fai * hi / ki;
    scat_h(i,3)=P;
%light holes to light holes,aborbtion---------------------------------
    ep = 1 / ( (1 / (eps_inf*eps_o)) - (1 / (eps_s*eps_o)) );
    N0 = 1 / (exp( (q*hw0)/(kb*T) ) - 1);
    c_const =(q/h)^3 * (hw0 * ml) / (4 * pi * ep);
    ab_const = c_const * N0;
    ei=delt_Ek*i;
    ef=ei+hw0;
    ki=sqrt(2*ml*delt_Ek*i*q/(h*h));
    kf=sqrt(2*ml*ef*q/(h*h));
    %kf=sqrt(ml/mh)*(ki^2+2*mh*w0/h);
    fai = log((ki + kf) / (max(ki, kf) - min(ki, kf)));        
    theta = (ki * ki + kf * kf) / (2 * ki * kf); 
    hi = (1 + 3 * theta * (theta - 1/fai)) / 4; 
    P=ab_const * fai * hi / ki;
    scat_l(i,3)=P;
%heavy holes to light holes,aborbtion---------------------------------
    ep = 1 / ( (1 / (eps_inf*eps_o)) - (1 / (eps_s*eps_o)) );
    N0 = 1 / (exp( (q*hw0)/(kb*T) ) - 1);
    c_const =(q/h)^3 * (hw0 * ml) / (4 * pi * ep);
    ab_const = c_const * N0;
    ei=delt_Ek*i;
    ef=ei+hw0;
    ki=sqrt(2*mh*delt_Ek*i*q/(h*h));
    kf=sqrt(2*ml*ef*q/(h*h));
    %kf=sqrt(ml/mh)*(ki^2+2*mh*w0/h);
    fai = log((ki + kf) / (max(ki, kf) - min(ki, kf)));        
    theta = (ki * ki + kf * kf) / (2 * ki * kf); 
    hi = 3 * (1 - theta * (theta - 1/fai)) / 4;
    P=ab_const * fai * hi / ki;
    scat_h(i,4)=P;
%light holes to heavy holes,aborbtion---------------------------------
    ep = 1 / ( (1 / (eps_inf*eps_o)) - (1 / (eps_s*eps_o)) );
    N0 = 1 / (exp( (q*hw0)/(kb*T) ) - 1);
    c_const =(q/h)^3 * (hw0 * mh) / (4 * pi * ep);
    ab_const = c_const * N0;
    ei=delt_Ek*i;
    ef=ei+hw0;
    ki=sqrt(2*ml*delt_Ek*i*q/(h*h));
    kf=sqrt(2*mh*ef*q/(h*h));
    %kf=sqrt(ml/mh)*(ki^2+2*mh*w0/h);
    fai = log((ki + kf) / (max(ki, kf) - min(ki, kf)));        
    theta = (ki * ki + kf * kf) / (2 * ki * kf); 
    hi = 3 * (1 - theta * (theta - 1/fai)) / 4;
    P=ab_const * fai * hi / ki;
    scat_l(i,4)=P;
%heavy holes to heavy holes,emission---------------------------------
    ep = 1 / ( (1 / (eps_inf*eps_o)) - (1 / (eps_s*eps_o)) );
    N0 = 1 / (exp( (q*hw0)/(kb*T) ) - 1);
    c_const =(q/h)^3 * (hw0 * mh) / (4 * pi * ep);
    em_const = c_const * (N0+1);
    ei=delt_Ek*i;
    ef=ei-hw0;
    if ef>de
        ki=sqrt(2*mh*delt_Ek*i*q/(h*h));
        kf=sqrt(2*mh*ef*q/(h*h));
        %kf=sqrt(ml/mh)*(ki^2+2*mh*w0/h);
        fai = log((ki + kf) / (max(ki, kf) - min(ki, kf)));        
        theta = (ki * ki + kf * kf) / (2 * ki * kf); 
        hi = (1 + 3 * theta * (theta - 1/fai)) / 4; 
        P=em_const * fai * hi / ki;
        scat_h(i,5)=P;
    else
        scat_h(i,5)=0;
    end
%light holes to light holes,emission---------------------------------
    ep = 1 / ( (1 / (eps_inf*eps_o)) - (1 / (eps_s*eps_o)) );
    N0 = 1 / (exp( (q*hw0)/(kb*T) ) - 1);
    c_const =(q/h)^3 * (hw0 * ml) / (4 * pi * ep);
    em_const = c_const * (N0+1);
    ei=delt_Ek*i;
    ef=ei-hw0;
    if ef>de
        ki=sqrt(2*ml*delt_Ek*i*q/(h*h));
        kf=sqrt(2*ml*ef*q/(h*h));
        %kf=sqrt(ml/mh)*(ki^2+2*mh*w0/h);
        fai = log((ki + kf) / (max(ki, kf) - min(ki, kf)));        
        theta = (ki * ki + kf * kf) / (2 * ki * kf); 
        hi = (1 + 3 * theta * (theta - 1/fai)) / 4; 
        P=em_const * fai * hi / ki;
        scat_l(i,5)=P;
    else
        scat_l(i,5)=0;
    end
%heavy holes to light holes,emission---------------------------------
    ep = 1 / ( (1 / (eps_inf*eps_o)) - (1 / (eps_s*eps_o)) );
    N0 = 1 / (exp( (q*hw0)/(kb*T) ) - 1);
    c_const =(q/h)^3 * (hw0 * ml) / (4 * pi * ep);
    em_const = c_const * (N0+1);
    ei=delt_Ek*i;
    ef=ei-hw0;
    if ef>de
        ki=sqrt(2*mh*delt_Ek*i*q/(h*h));
        kf=sqrt(2*ml*ef*q/(h*h));
        %kf=sqrt(ml/mh)*(ki^2+2*mh*w0/h);
        fai = log((ki + kf) / (max(ki, kf) - min(ki, kf)));        
        theta = (ki * ki + kf * kf) / (2 * ki * kf); 
        hi = 3 * (1 - theta * (theta - 1/fai)) / 4;
        P=em_const * fai * hi / ki;
        scat_h(i,6)=P;
    else
        scat_h(i,6)=0;
    end
%light holes to heavy holes,emission---------------------------------
    ep = 1 / ( (1 / (eps_inf*eps_o)) - (1 / (eps_s*eps_o)) );
    N0 = 1 / (exp( (q*hw0)/(kb*T) ) - 1);
    c_const =(q/h)^3 * (hw0 * mh) / (4 * pi * ep);
    em_const = c_const * (N0+1);
    ei=delt_Ek*i;
    ef=ei-hw0;
    if ef>de
        ki=sqrt(2*ml*delt_Ek*i*q/(h*h));
        kf=sqrt(2*mh*ef*q/(h*h));
        %kf=sqrt(ml/mh)*(ki^2+2*mh*w0/h);
        fai = log((ki + kf) / (max(ki, kf) - min(ki, kf)));        
        theta = (ki * ki + kf * kf) / (2 * ki * kf); 
        hi = 3 * (1 - theta * (theta - 1/fai)) / 4;
        P=em_const * fai * hi / ki;
        scat_l(i,6)=P;
    else
        scat_l(i,6)=0;
    end
end


totScath=scat_h(:,1)+scat_h(:,2)+scat_h(:,3)+scat_h(:,4)+scat_h(:,5)+scat_h(:,6)+scat_h(:,7)+scat_h(:,8)+scat_h(:,9)+scat_h(:,10)+scat_h(:,11)+scat_h(:,12);
totScatl=scat_l(:,1)+scat_l(:,2)+scat_l(:,3)+scat_l(:,4)+scat_l(:,5)+scat_l(:,6)+scat_l(:,7)+scat_l(:,8)+scat_l(:,9)+scat_l(:,10)+scat_l(:,11)+scat_l(:,12);

Gmh=max(totScath(:,1));
Gml=max(totScatl(:,1));

scatGaAs_hole(10,Ek_pts,2)=0;
scatGaAs_hole(1,:,1)=scat_h(:,1).';                   %1. Impurity,h-l
scatGaAs_hole(2,:,1)=scatGaAs_hole(1,:,1)+scat_h(:,2).';   %2. Impurity,h-h
scatGaAs_hole(3,:,1)=scatGaAs_hole(2,:,1)+scat_h(:,3).';   %3. Polar Optical Emission, h-h,a
scatGaAs_hole(4,:,1)=scatGaAs_hole(3,:,1)+scat_h(:,4).';   %4. Polar Optical Emission, h-l,a
scatGaAs_hole(5,:,1)=scatGaAs_hole(4,:,1)+scat_h(:,5).';   %5. Polar Optical Emission, h-h,3
scatGaAs_hole(6,:,1)=scatGaAs_hole(5,:,1)+scat_h(:,6).';   %6. Polar Optical Emission, h-l,e
scatGaAs_hole(7,:,1)=scatGaAs_hole(6,:,1)+scat_h(:,7).';   %7. NPOP, h-h,a
scatGaAs_hole(8,:,1)=scatGaAs_hole(7,:,1)+scat_h(:,8).';   %8. NPOP, h-h,e
scatGaAs_hole(9,:,1)=scatGaAs_hole(8,:,1)+scat_h(:,9).';   %9. NPOP, h-l,a
scatGaAs_hole(10,:,1)=scatGaAs_hole(9,:,1)+scat_h(:,10).';   %10. NPOP, h-l,e
scatGaAs_hole(11,:,1)=scatGaAs_hole(10,:,1)+scat_h(:,11).';   %10. NPOP, h-l,e
scatGaAs_hole(12,:,1)=scatGaAs_hole(11,:,1)+scat_h(:,12).';
scatGaAs_hole(:,:,1)=scatGaAs_hole(:,:,1)./Gmh;

scatGaAs_hole(1,:,2)=scat_l(:,2).';                   %1. Impurity,h-l
scatGaAs_hole(2,:,2)=scatGaAs_hole(1,:,2)+scat_l(:,2).';   %2. Impurity,h-h
scatGaAs_hole(3,:,2)=scatGaAs_hole(2,:,2)+scat_l(:,3).';   %3. Polar Optical Emission, h-h,a
scatGaAs_hole(4,:,2)=scatGaAs_hole(3,:,2)+scat_l(:,4).';   %4. Polar Optical Emission, h-l,a
scatGaAs_hole(5,:,2)=scatGaAs_hole(4,:,2)+scat_l(:,5).';   %5. Polar Optical Emission, h-h,3
scatGaAs_hole(6,:,2)=scatGaAs_hole(5,:,2)+scat_l(:,6).';   %6. Polar Optical Emission, h-l,e
scatGaAs_hole(7,:,2)=scatGaAs_hole(6,:,2)+scat_l(:,7).';   %7. NPOP, h-h,a
scatGaAs_hole(8,:,2)=scatGaAs_hole(7,:,2)+scat_l(:,8).';   %8. NPOP, h-h,e
scatGaAs_hole(9,:,2)=scatGaAs_hole(8,:,2)+scat_l(:,9).';   %9. NPOP, h-l,a
scatGaAs_hole(10,:,2)=scatGaAs_hole(9,:,2)+scat_l(:,10).';   %10. NPOP, h-l,e
scatGaAs_hole(11,:,2)=scatGaAs_hole(10,:,2)+scat_l(:,11).';
scatGaAs_hole(12,:,2)=scatGaAs_hole(11,:,2)+scat_l(:,12).';
scatGaAs_hole(:,:,2)=scatGaAs_hole(:,:,2)./Gml;

