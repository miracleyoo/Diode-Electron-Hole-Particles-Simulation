clear all;
clc;

%-------Simulation Input Settings------------------------------------------
tsteps=20000;
dt=10e-15;
%include hole
dope_type=1;
cdop=1e23;
cdop2=2e20;
cdop3=0;%for substrate0.1
hd=3;%highly doped mear contact
max_particles=40e4;%55e4;
de=0.002;
T=300;
Vp_all=[0];
Vn_all=0;
POISSON_ITER_MAX=12000;
ppc=8;

path_str='C:\Users\sunz1\Desktop\boeing_dropbox\nick_matlab_code\hole_scattering\pn_junction_data\1';

%--------Simulation Settings-----------------------------------------------
dx=1e-8;
dy=1e-8;

%--------Device Geomrety---------------------------------------------------
Ttot=0.08e-6;
Ltot=0.4e-6;
Lp=Ltot/2;
Ln=Ltot/2;

%---------Constants--------------------------------------------------------
bk=1.38066e-23;                 %Boltzmann's Constant 
q=1.60219e-19;                  %Charge of Electron
h=1.05459e-34;                  %Planck's Constant (/2pi)
emR=9.10953e-31;                %Mass of Resting Electron
eps_o=8.85419e-12;              %Vacuum Permittivity

%---------GaAs Specific Constants------------------------------------------
Eg=1.424;                       %Band Gap for GaAs
Egg=0;                          %Energy difference between two generate Gamma Bands
Egl=0.29;                       %Energy difference between Gamma and L Bands
eC=[Egg Egl];
emG=0.067*emR;                  %Mass of Electron in Gamma Band
emL=0.350*emR;                  %Mass of Electron in L Band
emh=0.62*emR;
eml=0.074*emR;
eM=[emG,emL,emh,eml];
alpha_G=(1/Eg)*(1-emG/emR)^2;   %Alpha for Non-Parabolicity of Gamma Band
alpha_L=(1/(Eg+Egl))*(1-emL/emR)^2;%Alpha for Non-Parabolicty of L Band
alpha=[alpha_G,alpha_L,0,0];%for holes, alpha=0(assume parabolic)

eps_stat=12.9*eps_o;            %Static Permittivity for GaAs
eps_inf=10.92*eps_o;            %Optical Permittivity for GaAs
eps_p=1/((1/eps_inf)-(1/eps_stat)); 

qD=sqrt(q*q*cdop/(eps_stat*bk*T)); %Inverse Debye Length

ni=1.8e12;%GaAs at 300k
contact_potential=bk*T/q*log(cdop*cdop/ni^2);

hw0=0.03536;
hwij=0.03;
hwe=hwij;

%inverse band mass parameters
A=-7.65;
B=-4.82;
C=7.7;
g100=B/A;
g111=sqrt((B/A)^2+(C/A)^2/3);

%---------Create Scattering Table------------------------------------------
[scatGaAs,GmG,GmL]=make_GaAs_scatTable(T,0,de,2,cdop);%change 2 to Vmax, 108
[scatGaAs_hole,Gmh,Gml]=make_GaAs_hole_scatTable_v2(T,de,2,cdop);
Gm=[GmG,GmL,Gmh,Gml];
Gm_max=max(Gmh,Gml);

%------------related to configuration--------------------------------------
nx1=round(Ltot/dx)+1;
nx=nx1-1;

ny1=round(Ttot/dy)+1;
ny=ny1-1;

bottom_pts=2:nx;
left_pts=1:nx1:nx1*ny1;
right_pts=nx1:nx1:nx1*ny1;
top_pts=nx1*ny1-nx1+2:nx1*ny1-1;

p_icpg(length(left_pts))=0;
n_icpg(length(right_pts))=0;

% tic
for a=1:length(Vp_all)
    for b=1:length(Vn_all)
        Vp=Vp_all(a)-contact_potential/2;
        Vn=Vn_all(b)+contact_potential/2;
        
        %---------Initialize Particles---------------------------------------------
        [particles,valley,bg_charge,cpsp,N,p_icpg,n_icpg,xmax,ymax]=pn_init_v2(max_particles,dope_type,dx,dy,Ltot,nx1,ny1,cdop,ppc,bk,T,q,h,alpha,eM,Gm,Lp,A,B,C,emR,hd);
        %grids start at 1, and position starts at 0

        %--------Initial Charge/Field Computations---------------------------------
        [charge_p,charge_n]=pn_charge_v2(particles,valley,nx1,ny1,dx,dy,max_particles,cpsp);%¾»µç×ÓÊý
        phi=zeros(ny1,nx1);
        [fx,fy,phi,k]=pn_poisson_v5(dx,dy,nx1,ny1,eps_stat,q,charge_p,charge_n,bg_charge,phi,Vp,Vn);

        number=zeros(tsteps,4);
        
        %calculate current
        p_enter=zeros(tsteps,1);
        p_exit=zeros(tsteps,1);
        p_real=zeros(tsteps,1);
        n_enter=zeros(tsteps,1);
        n_exit=zeros(tsteps,1);
        n_real=zeros(tsteps,1);
        current_n=zeros(tsteps,1);
        current_p=zeros(tsteps,1);

        %record energy of every particle at every time step
      
        t_poisson=0;
        input=zeros(ny1,nx1,2,tsteps);
        ii=1;
        
        %------------Main Loop----------------------------------------------------
        for ti=1:tsteps
            tic;
            p_temp=0;
            n_temp=0;
            t=(ti-1)*dt;
            %EMC Routine
            tdt=t+dt;
            for n=1:max_particles
                if valley(n,1)~=9
                    ts=particles(n,4);
                    t1=t;
                    while ts<tdt
                        tau=ts-t1;
                        %Drift------------------------------
                        iv=valley(n,1);
                        if iv~=9
                            if iv==1||iv==2
                                kx=particles(n,1);
                                ky=particles(n,2);
                                kz=particles(n,3);
                                x=particles(n,5);
                                y=particles(n,6);

                                i=min(floor(y/dy)+1,ny1);
                                j=min(floor(x/dx)+1,nx1);
                                i=max(i,1);
                                j=max(j,1);%928

                                dkx=-(q/h)*fx(i,j)*tau;%the field can be more accurate by interprolate
                                dky=-(q/h)*fy(i,j)*tau;
                                sk=kx*kx+ky*ky+kz*kz;
                                gk=(h*h/(2*eM(iv)))*sk*(1/q);

                                x=x+(h/eM(iv))*(kx+0.5*dkx)/(sqrt(1+4*alpha(iv)*gk))*tau;
                                y=y+(h/eM(iv))*(ky+0.5*dky)/(sqrt(1+4*alpha(iv)*gk))*tau;

                                particles(n,1)=kx+dkx;  %kx
                                particles(n,2)=ky+dky;  %ky
                            else
                                kx=particles(n,1);
                                ky=particles(n,2);
                                kz=particles(n,3);
                                x=particles(n,5);
                                y=particles(n,6);

                                i=min(floor(y/dy)+1,ny1);
                                j=min(floor(x/dx)+1,nx1);
                                i=max(i,1);
                                j=max(j,1);

                                dkx=(q/h)*fx(i,j)*tau;
                                dky=(q/h)*fy(i,j)*tau;                                            

                                if iv==3%heavy                                    
                                    kf=sqrt((kx+dkx)^2+(ky+dky)^2+kz*kz);        
                                    cos_theta=kz/kf;
                                    sin_theta=sqrt(1-cos_theta^2);
                                    sin_phi=ky/kf/sin_theta;
                                    cos_phi=(kx+dkx)/kf/sin_theta;
                                    g=((B/A)^2+(C/A)^2*(sin_theta^2*cos_theta^2+sin_theta^4*cos_phi^2*sin_phi^2))^0.5;   
                                    mh=emR/(abs(A)*(1-g));
                                    x=x+(h/mh)*(kx+0.5*dkx)*tau;
                                    y=y+(h/mh)*(ky)*tau;
                                    particles(n,1)=kx+dkx;  %kx
                                    particles(n,2)=ky+dky;  %ky
                                elseif iv==4
                                    kf=sqrt((kx+dkx)^2+(ky+dky)^2+kz*kz);
                                    cos_theta=kz/kf;
                                    sin_theta=sqrt(1-cos_theta^2);
                                    sin_phi=ky/kf/sin_theta;
                                    cos_phi=(kx+dkx)/kf/sin_theta;
                                    g=((B/A)^2+(C/A)^2*(sin_theta^2*cos_theta^2+sin_theta^4*cos_phi^2*sin_phi^2))^0.5;   
                                    ml=emR/(abs(A)*(1+g));
                                    ef=(h*h*abs(A)/(2*emR))*kf^2*(1/q)*(1+g);
                                    x=x+(h/ml)*(kx+0.5*dkx)*tau;
                                    y=y+(h/ml)*(ky)*tau;
                                    particles(n,1)=kx+dkx;  %kx
                                    particles(n,2)=ky+dky;  %ky
                                end
                            end   
                            %Boundary Condition-----the former change is incorrect, only kx or ky one has to change 921----------------
                            if x < 0
                                valley(n,1)=9;
                                if iv==1||iv==2
                                    p_temp=p_temp-1;%p_temp count how mant positive particles are deleted when drifting
                                else
                                    p_temp=p_temp+1;
                                end
                            elseif x > xmax
                                valley(n,1)=9; 
                                if iv==1||iv==2
                                    n_temp=n_temp+1;%n_temp count how mant negative particles are deleted when drifting
                                else
                                    n_temp=n_temp-1;
                                end
                            end

                            if y > ymax
                                y=ymax-(y-ymax);
                                particles(n,2)=-particles(n,2);
                            elseif y < 0
                                y=-y;
                                particles(n,2)=-particles(n,2);
                            end

                            particles(n,5)=x;       %x
                            particles(n,6)=y;       %y

                            %Scatter----------------------------
                            if valley(n,1)~=9
                                [particle,valley(n,1)]=pn_scat_v2(particles(n,:),valley(n,1),scatGaAs,scatGaAs_hole,de,q,h,eM,alpha,qD,hw0,A,B,C,emR,n,hwij,Egl,Egg,hwe,g100,g111);
                                particles(n,:)=particle(1,:);
                            end

                            t1=ts;
                            ts=t1-log(rand())/Gm(iv);%can be more accurate
                        else
                            ts=tdt;
                        end
                    end

                    tau=tdt-t1;
                    %Drift------------------------------------
                    iv=valley(n,1);
                    if iv~=9                        
%                         tic
                        if iv==1||iv==2
                            kx=particles(n,1);
                            ky=particles(n,2);
                            kz=particles(n,3);
                            x=particles(n,5);
                            y=particles(n,6);

                            i=min(floor(y/dy)+1,ny1);
                            j=min(floor(x/dx)+1,nx1);
                            i=max(i,1);
                            j=max(j,1);%928

                            dkx=-(q/h)*fx(i,j)*tau;%the field can be more accurate by interprolate
                            dky=-(q/h)*fy(i,j)*tau;
                            sk=kx*kx+ky*ky+kz*kz;
                            gk=(h*h/(2*eM(iv)))*sk*(1/q);

                            x=x+(h/eM(iv))*(kx+0.5*dkx)/(sqrt(1+4*alpha(iv)*gk))*tau;
                            y=y+(h/eM(iv))*(ky+0.5*dky)/(sqrt(1+4*alpha(iv)*gk))*tau;

                            particles(n,1)=kx+dkx;  %kx
                            particles(n,2)=ky+dky;  %ky
                        else
                            t0=cputime;
                            kx=particles(n,1);
                            ky=particles(n,2);
                            kz=particles(n,3);
                            x=particles(n,5);
                            y=particles(n,6);

                            i=min(floor(y/dy)+1,ny1);
                            j=min(floor(x/dx)+1,nx1);
                            i=max(i,1);
                            j=max(j,1);

                            dkx=(q/h)*fx(i,j)*tau;
                            dky=(q/h)*fy(i,j)*tau;     

                            if iv==3%heavy                                    
                                kf=sqrt((kx+dkx)^2+(ky+dky)^2+kz*kz);        
                                cos_theta=kz/kf;
                                sin_theta=sqrt(1-cos_theta^2);
                                sin_phi=ky/kf/sin_theta;
                                cos_phi=(kx+dkx)/kf/sin_theta;
                                g=((B/A)^2+(C/A)^2*(sin_theta^2*cos_theta^2+sin_theta^4*cos_phi^2*sin_phi^2))^0.5;   
                                mh=emR/(abs(A)*(1-g));
                                x=x+(h/mh)*(kx+0.5*dkx)*tau;
                                y=y+(h/mh)*(ky)*tau;
                                particles(n,1)=kx+dkx;  %kx
                                particles(n,2)=ky+dky;  %ky
                            elseif iv==4
                                kf=sqrt((kx+dkx)^2+(ky+dky)^2+kz*kz);
                                cos_theta=kz/kf;
                                sin_theta=sqrt(1-cos_theta^2);
                                sin_phi=ky/kf/sin_theta;
                                cos_phi=(kx+dkx)/kf/sin_theta;
                                g=((B/A)^2+(C/A)^2*(sin_theta^2*cos_theta^2+sin_theta^4*cos_phi^2*sin_phi^2))^0.5;   
                                ml=emR/(abs(A)*(1+g));
                                ef=(h*h*abs(A)/(2*emR))*kf^2*(1/q)*(1+g);
                                x=x+(h/ml)*(kx+0.5*dkx)*tau;
                                y=y+(h/ml)*(ky)*tau;
                                particles(n,1)=kx+dkx;  %kx
                                particles(n,2)=ky+dky;  %ky
                            end

                        end   
                        %Boundary Condition-----the former change is incorrect, only kx or ky one has to change 921----------------
                        if x < 0
                            valley(n,1)=9;
                            if iv==1||iv==2
                                p_temp=p_temp-1;%p_temp count how mant positive particles are deleted when drifting
                            else
                                p_temp=p_temp+1;
                            end
                        elseif x > xmax
                            valley(n,1)=9; 
                            if iv==1||iv==2
                                n_temp=n_temp+1;%n_temp count how mant negative particles are deleted when drifting
                            else
                                n_temp=n_temp-1;
                            end
                        end

                        if y > ymax
                            y=ymax-(y-ymax);
                            particles(n,2)=-particles(n,2);
                        elseif y < 0
                            y=-y;
                            particles(n,2)=-particles(n,2);
                        end

                        particles(n,5)=x;       %x
                        particles(n,6)=y;       %y 

                        particles(n,4)=ts;      %ts
                    end
                end
            end              

            %Renew--------------------
            [particles,valley,p_added,n_added,number]=pn_renew_v6(particles,valley,Ttot,dx,dy,nx1,ny1,max_particles,p_icpg,n_icpg,bk,T,q,h,alpha,eM,emR,Gm,tdt,left_pts,right_pts,Ltot,A,B,C,ti,number,hd);
            p_real(ti,1)=p_added-p_temp;%p_added is how many positive particles are injected in renew
            n_real(ti,1)=n_added-n_temp;%n_added is how many negative particles are injected in renew
            %p_real is how many positive particles are injected in ti
            %n_real is how many negative particles are injected in ti

            %Charge Computation-------
            [charge_p,charge_n]=pn_charge_v2(particles,valley,nx1,ny1,dx,dy,max_particles,cpsp);

            %Poisson Calculation
            tic
            [fx,fy,phi,k]=pn_poisson_v5(dx,dy,nx1,ny1,eps_stat,q,charge_p,charge_n,bg_charge,phi,Vp,Vn);
            t_poisson=toc+t_poisson;

            %net cathode current----
            if ti==1
                current_n(ti,1)=current_n(ti,1);
            else
                current_n(ti,1)=current_n(ti-1,1)+n_real(ti,1);
            end
            %net anode current----
            if ti==1
                current_p(ti,1)=current_p(ti,1);
            else
                current_p(ti,1)=current_p(ti-1,1)+p_real(ti,1);
            end

            index_v1=valley(:,1)==1;
            index_v1=index_v1(:,1).*(1:max_particles).';
            index_v1=index_v1(index_v1~=0);
% 
            index_v2=valley(:,1)==2;
            index_v2=index_v2(:,1).*(1:max_particles).';
            index_v2=index_v2(index_v2~=0);
            
            index_v3=valley(:,1)==3;
            index_v3=index_v3(:,1).*(1:max_particles).';
            index_v3=index_v3(index_v3~=0);
            
            index_v4=valley(:,1)==4;
            index_v4=index_v4(:,1).*(1:max_particles).';
            index_v4=index_v4(index_v4~=0);
            
            particle_num(ti,:)=[length(index_v1),length(index_v2),length(index_v3),length(index_v4)];
            jj='progress';%914
            fprintf('%s %f\n',jj,ti);%914
        end
        
        aa=Vp_all(a)*100;
        bb=Vn_all(b)*100;
        save(['Vp=',num2str(aa),'Vn=',num2str(bb),'dx=',num2str(dx*1e9),'nm']);%bceause it can't save files like 0.1 
    end
end