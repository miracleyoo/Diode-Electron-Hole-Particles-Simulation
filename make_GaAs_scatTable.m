function [scatGaAs,GmG,GmL]=make_GaAs_scatTable(T,spherical_only,de,Vmax,cimp)

%--------Electron Energy Steps for Formulas/Graphs-------------------------
delt_Ek=de;                     %Electron Energy Step (default 0.002);
Ek_pts=Vmax/de;                 %Number of Sample Points for Electron Energy   
eV_axis=(1:Ek_pts)*delt_Ek;
scatG(Ek_pts,6)=0;              %Matrix holding scattering rates for Gamma Valley
scatL(Ek_pts,8)=0;              %Matrix holding scattering rates for L Valley

%---------General Constants------------------------------------------------
q=1.60219e-19;                  %Charge of Electron
bk=1.38066e-23;                 %Boltzmann's Constant
h=1.05459e-34;                  %Planck's Constant (/2pi)
emR=9.10953e-31;                %Mass of Electron at Rest
eps_o=8.85419e-12;              %Vacuum Permittivity

%---------GaAs Specific Constants------------------------------------------
Eg=1.424;                           %Band Gap for GaAs
Egg=0;                              %Energy difference between two generate Gamma Bands
Egl=0.29;                           %Energy difference between Gamma and L Bands
emG=0.067*emR;                      %Mass of Electron in Gamma Band
emL=0.350*emR;                      %Mass of Electron in L Band

kconst_G=sqrt(2*emG)/h;             %Constant part of k for Parabolic Gamma Band
kconst_L=sqrt(2*emL)/h;             %Constant part of k for Parabolic L Band

alpha_G=(1/Eg)*(1-emG/emR)^2;       %Alpha for Non-Parabolicity of Gamma Band
alpha_L=(1/(Eg+Egl))*(1-emL/emR)^2; %Alpha for Non-Parabolicty of L Band

eps_stat=12.9*eps_o;                %Static Permittivity for GaAs
eps_inf=10.92*eps_o;                %Optical Permittivity for GaAs
eps_p=1/((1/eps_inf)-(1/eps_stat));

%---------GaAs Constants for Acoustic Phonon Scattering--------------------
rho=5360;                           %Crystal Density
sv=5240;                            %Sound Velocity of Longitudinal Elastic Waves
cl=rho*sv*sv;                       %Elastic Constant of Material
adp=7*q;                            %Acoustic Deformation Potential
acoustic_const=(2*pi*adp*adp*bk*T)/(h*cl);  %Constant in Acoustic Scattering Formula

%--------GaAs Constants for Polar Optical Phonon Scattering----------------
hwo=0.03536;                    %In eV -> Need 0.03536*q for calculations
hwoq=hwo*q;
wo=hwoq/h;                      %Transverse Polar Optical Frequency??
no=1/(exp((hwoq)/(bk*T))-1);    %Number of Phonons with Mode o at Temp T
pop_constG=((q*q*wo)/(8*pi*eps_p))*kconst_G;%Constant in Polar Optical Scattering Formula
pop_constL=((q*q*wo)/(8*pi*eps_p))*kconst_L;%Constant in Polar Optical Scattering Formula

%-------GaAs Constants for Intervalley Scattering - Non Polar Optical Phonon
hwij=0.03;                      %In eV 
hwijq=hwij*q;
wij=hwijq/h;
dij=1e11*q;                     %Intervalley Deformation Potential
Zgl=4;                          %Number of Equivalent Valleys for Scattering from Gamma to L Band
Zlg=1;                          %Number of Equivalent Valleys for Scattering from L Band to Gamma Band
nij=1/(exp(hwijq/(bk*T))-1);    %Number of Phonons with Mode ij at Temp T
npop_constGL=(pi*dij*dij*Zgl)/(rho*wij);
npop_constLG=(pi*dij*dij*Zlg)/(rho*wij);

%------GaAs Constants for Intravalley Scattering - Non Polar Optical Phonon
hwll=hwij;
hwllq=hwll*q;
Zll=Zgl-1;
dll=dij;
wll=hwllq/h;
nll=1/(exp(hwllq/(bk*T))-1);
npop_constLL=(pi*dll*dll*Zll)/(rho*wll);

%-------GaAs Constants for Impurity Scattering-----------------------------
qD=sqrt(q*q*cimp/(eps_stat*bk*T)); %Inverse Debye Length
NI=cimp;
imp_const=(2*pi*NI*q*q*q*q)/(h*eps_stat*eps_stat);


%------------Acoustic Phonon Scattering------------------------------------
for i=1:Ek_pts
    Ek=delt_Ek*i;
    
    %-------Gamma Valley---------------------
    %For Spherical Parabolic
    if spherical_only == 1
        N_Ek=(((2*emG)^(3/2))/(4*pi*pi*h*h*h))*sqrt(Ek)*sqrt(q);                
    %For Spherical + Non-Parabolicty Factor
    else
        gamma_part=sqrt(Ek*(1+alpha_G*Ek))*(1+2*alpha_G*Ek)*sqrt(q);    %Gamma Part of Density of States Formula
        N_Ek=(((2*emG)^(3/2))/(4*pi*pi*h*h*h))*gamma_part;              %Density of States Calculation For Spherical Non Parabolic
    end
    scatG(i,1)=acoustic_const*N_Ek;
    
    %--------L Valley------------------------
    if spherical_only==1
        N_Ek=(((2*emL)^(3/2))/(4*pi*pi*h*h*h))*sqrt(Ek)*sqrt(q);
    else
        gamma_part=sqrt(Ek*(1+alpha_L*Ek))*(1+2*alpha_L*Ek)*sqrt(q);    %Gamma Part of Density of States Formula
        N_Ek=(((2*emL)^(3/2))/(4*pi*pi*h*h*h))*gamma_part;              %Density of States Calculation For Spherical Non Parabolic
    end
    scatL(i,1)=acoustic_const*N_Ek;
end

%-----------Polar Optical Phonon Scattering--------------------------------
for i=1:Ek_pts
    Ek=delt_Ek*i;
    
    %Gamma Band---------------------------------------
    %CaseI: 1+hwo/Ek (Absorption)
    qmax=sqrt(Ek)+sqrt(Ek+hwo);
    qmin=sqrt(Ek+hwo)-sqrt(Ek);
    scatG(i,2)=((pop_constG*no*sqrt(Ek*q))/(Ek*q))*log(qmax/qmin);
    
    %CaseII: 1-hwo/Ek (Emission)
    if Ek-hwo > 0
        qmax=sqrt(Ek)+sqrt(Ek-hwo);
        qmin=sqrt(Ek)-sqrt(Ek-hwo);
        scatG(i,3)=(pop_constG*(no+1)*sqrt(Ek*q)/(Ek*q))*log(qmax/qmin);
    else
        scatG(i,3)=0;
    end
    
    %L Band------------------------------------------
    %CaseI: 1+hwo/Ek (Absorption)
    qmax=sqrt(Ek)+sqrt(Ek+hwo);
    qmin=sqrt(Ek+hwo)-sqrt(Ek);
    scatL(i,2)=((pop_constL*no*sqrt(Ek*q))/(Ek*q))*log(qmax/qmin);
    
    %CaseII: 1-hwo/Ek (Emission)
    if Ek-hwo > 0
        qmax=sqrt(Ek)+sqrt(Ek-hwo);
        qmin=sqrt(Ek)-sqrt(Ek-hwo);
        scatL(i,3)=(pop_constL*(no+1)*sqrt(Ek*q)/(Ek*q))*log(qmax/qmin);
    else
        scatL(i,3)=0;
    end
end

%---------Intervalley Scattering - Non Polar Optical Phonons---------------
for i=1:Ek_pts
    Ek=delt_Ek*i;
    %Gamma to L------------------------------------------------------------
    %CaseI (Absorption):
    Ek2=Ek+hwij-(Egl-Egg);
    if Ek2 > 0
        %For Spherical Parabolic
        if spherical_only == 1
            N_Ek2=(((2*emL)^(3/2))/(4*pi*pi*h*h*h))*sqrt(Ek2)*sqrt(q);
        %For Spherical Non-Parabolic
        else
            gamma_part=sqrt(Ek2*(1+alpha_L*Ek2))*(1+2*alpha_L*Ek2)*sqrt(q); %Gamma Part of Density of States Formula 
            N_Ek2=(((2*emL)^(3/2))/(4*pi*pi*h*h*h))*gamma_part;    %Density of States Calculation
        end
        scatG(i,4)=npop_constGL*N_Ek2*nij;
    else
        scatG(i,4)=0;
    end
    
    %CaseII (Emission): 
    Ek2=Ek-hwij-(Egl-Egg);
    if Ek2 > 0
        %For Spherical Parabolic
        if spherical_only == 1
            N_Ek2=(((2*emL)^(3/2))/(4*pi*pi*h*h*h))*sqrt(Ek2)*sqrt(q);
        %For Spherical Non-Parabolic
        else
            gamma_part=sqrt(Ek2*(1+alpha_L*Ek2))*(1+2*alpha_L*Ek2)*sqrt(q); %Gamma Part of Density of States Formula 
            N_Ek2=(((2*emL)^(3/2))/(4*pi*pi*h*h*h))*gamma_part;    %Density of States Calculation
        end
        scatG(i,5)=npop_constGL*N_Ek2*(nij+1);
    else
        scatG(i,5)=0;
    end
    
    %L to Gamma------------------------------------------------------------
    %CaseI (Absorption):
    Ek2=Ek+hwij-(Egg-Egl);
    if Ek2 > 0
        %For Spherical Parabolic
        if spherical_only == 1
            N_Ek2=(((2*emG)^(3/2))/(4*pi*pi*h*h*h))*sqrt(Ek2)*sqrt(q);
        else
            gamma_part=sqrt(Ek2*(1+alpha_G*Ek2))*(1+2*alpha_G*Ek2)*sqrt(q); %Gamma Part of Density of States Formula 
            N_Ek2=(((2*emG)^(3/2))/(4*pi*pi*h*h*h))*gamma_part;    %Density of States Calculation
        end
        scatL(i,4)=npop_constLG*N_Ek2*nij;
    else
        scatL(i,4)=0;
    end
    
    %CaseII (Emission): 
    Ek2=Ek-hwij-(Egg-Egl);
    if Ek2 > 0
        %For Spherical Parabolic
        if spherical_only == 1
            N_Ek2=(((2*emG)^(3/2))/(4*pi*pi*h*h*h))*sqrt(Ek2)*sqrt(q);
        else
            gamma_part=sqrt(Ek2*(1+alpha_G*Ek2))*(1+2*alpha_G*Ek2)*sqrt(q); %Gamma Part of Density of States Formula 
            N_Ek2=(((2*emG)^(3/2))/(4*pi*pi*h*h*h))*gamma_part;    %Density of States Calculation
        end
        scatL(i,5)=npop_constLG*N_Ek2*(nij+1);
    else
        scatL(i,5)=0;
    end
end

%---------Intravalley Scattering - Non Polar Optical Phonon----------------
%For L Valley Only
for i=1:Ek_pts
    Ek=delt_Ek*i;
    
    %CaseI
    Ek2=Ek+hwll;
    if spherical_only==1
        N_Ek2=(((2*emL)^(3/2))/(4*pi*pi*h*h*h))*sqrt(Ek2)*sqrt(q);
    else
        gamma_part=sqrt(Ek2*(1+alpha_L*Ek2))*(1+2*alpha_L*Ek2)*sqrt(q);
        N_Ek2=(((2*emL)^(3/2))/(4*pi*pi*h*h*h))*gamma_part;
    end
    scatL(i,6)=npop_constLL*N_Ek2*nll;
    
    %CaseII
    Ek2=Ek-hwll;
    if Ek2 > 0
        if spherical_only==1
            N_Ek2=(((2*emL)^(3/2))/(4*pi*pi*h*h*h))*sqrt(Ek2)*sqrt(q);
        else
            gamma_part=sqrt(Ek2*(1+alpha_L*Ek2))*(1+2*alpha_L*Ek2)*sqrt(q);
            N_Ek2=(((2*emL)^(3/2))/(4*pi*pi*h*h*h))*gamma_part;
        end
        scatL(i,7)=npop_constLL*N_Ek2*(1+nll);
    else
        scatL(i,7)=0;
    end
end

%--------Impurity Scattering-----------------------------------------------
for i=1:Ek_pts
    %Gamma Band---------------------------------
    Ek=delt_Ek*i;
    if spherical_only==1
        kpart=sqrt(2*emG*Ek*q)/h;
        denom=qD*qD*(4*kpart*kpart+qD*qD);
        N_Ek=(((2*emG)^(3/2))/(4*pi*pi*h*h*h))*sqrt(Ek)*sqrt(q);
    else
        kpart=sqrt(Ek*(1+alpha_G*Ek))*sqrt(q);
        denom=qD*qD*(4*kconst_G*kconst_G*kpart*kpart+qD*qD);
        gamma_part=sqrt(Ek*(1+alpha_G*Ek))*(1+2*alpha_G*Ek)*sqrt(q);
        N_Ek=(((2*emG)^(3/2))/(4*pi*pi*h*h*h))*gamma_part;
    end
    scatG(i,6)=imp_const*(1/denom)*N_Ek;
    
    %L Band-------------------------------------
    Ek=delt_Ek*i;
    if spherical_only==1
        kpart=sqrt(2*emL*Ek*q)/h;
        denom=qD*qD*(4*kpart*kpart+qD*qD);
        N_Ek=(((2*emL)^(3/2))/(4*pi*pi*h*h*h))*sqrt(Ek)*sqrt(q);
    else
        kpart=sqrt(Ek*(1+alpha_L*Ek))*sqrt(q);
        denom=qD*qD*(4*kconst_L*kconst_L*kpart*kpart+qD*qD);
        gamma_part=sqrt(Ek*(1+alpha_L*Ek))*(1+2*alpha_L*Ek)*sqrt(q);
        N_Ek=(((2*emL)^(3/2))/(4*pi*pi*h*h*h))*gamma_part;
    end
    scatL(i,8)=imp_const*(1/denom)*N_Ek;
end


%Total Scattering Rate = Sum(All Scatter Mechanisms)
totScatG=scatG(:,1)+scatG(:,2)+scatG(:,3)+scatG(:,4)+scatG(:,5)+scatG(:,6);
totScatL=scatL(:,1)+scatL(:,2)+scatL(:,3)+scatL(:,4)+scatL(:,5)+scatL(:,6)+scatL(:,7)+scatL(:,8);


GmG=max(totScatG(:,1));
GmL=max(totScatL(:,1));

scatGaAs(8,Ek_pts,2)=0;
scatGaAs(1,:,1)=scatG(:,1).';                   %1. Acoustic
scatGaAs(2,:,1)=scatGaAs(1,:,1)+scatG(:,2).';   %2. Polar Optical Absorption +hwo
scatGaAs(3,:,1)=scatGaAs(2,:,1)+scatG(:,3).';   %3. Polar Optical Emission -hwo
scatGaAs(4,:,1)=scatGaAs(3,:,1)+scatG(:,4).';   %4. NPOP (Intervalley) Absorption +hwij
scatGaAs(5,:,1)=scatGaAs(4,:,1)+scatG(:,5).';   %5. NPOP (Intervalley) Emission -hwij
scatGaAs(6,:,1)=scatGaAs(5,:,1)+scatG(:,6).';   %6. Impurity
scatGaAs(:,:,1)=scatGaAs(:,:,1)./GmG;

scatGaAs(1,:,2)=scatL(:,1).';                   %1. Acoustic
scatGaAs(2,:,2)=scatGaAs(1,:,2)+scatL(:,4).';   %2. NPOP (Intervalley) Absorption +hwij
scatGaAs(3,:,2)=scatGaAs(2,:,2)+scatL(:,5).';   %3. NPOP (Intervalley) Emission -hwij
scatGaAs(4,:,2)=scatGaAs(3,:,2)+scatL(:,2).';   %4. POP Absorption +hwo
scatGaAs(5,:,2)=scatGaAs(4,:,2)+scatL(:,3).';   %5. POP Emission -hwo
scatGaAs(6,:,2)=scatGaAs(5,:,2)+scatL(:,6).';   %6. NPOP (Intravalley) Absorption +hwe
scatGaAs(7,:,2)=scatGaAs(6,:,2)+scatL(:,7).';   %7. NPOP (Intravalley) Emission -hwe
scatGaAs(8,:,2)=scatGaAs(7,:,2)+scatL(:,8).';   %8. Impurity 
scatGaAs(:,:,2)=scatGaAs(:,:,2)./GmL;