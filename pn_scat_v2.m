function [particle,valleyf]=pn_scat_v2(particle,valley,scatGaAs,scatGaAs_hole,de,q,h,eM,alpha,qD,hw0,A,B,C,m0,n,hwij,Egl,Egg,hwe,g100,g111)
% function [particle,valleyf]=bulk_scat(particle,valley,scatTable,eM,alpha,xs,xd,qD,hwo,hwij,hwe,Egl,Egg)
%particle: vector for particle properties (one row of particles matrix)
%valley: valley index for particle (1:Gamma, 2:L, 9:Delete)
%scatTable: table with precalculated scattering rates
%eM: effective electron mass vector
%alpha: non-parabolicity factor (vector for Gamma and L)
%xs: x-coord of end of cathode region
%xd: x-coord of start of anode region
%qD: inverse debye length

kx=particle(1,1);
ky=particle(1,2);
kz=particle(1,3);
a=h*h*abs(A)/(2*m0);
b=qD;

skx=kx*kx;
sky=ky*ky;
skz=kz*kz;
sk=abs(skx+sky+skz);
ki=sqrt(sk);

iemax=length(scatGaAs(1,:,1));%how many Ek_pts=Vmax/de;

kxf=kx;
kyf=ky;
kzf=kz;
valleyf=valley;

if sk ~= 0  
    gk=(h*h/(2*eM(valley)))*sk*(1/q);
    if valley==1||valley==2
        ei=(sqrt(1+4*alpha(valley)*gk)-1)/(2*alpha(valley));
        particle(1,7)=ei;
    elseif valley==3
        cos_theta=kz/ki;
        sin_theta=sqrt(1-cos_theta^2);
        cos_phi=kx/(ki*sin_theta);
        sin_phi=ky/(ki*sin_theta);
        g1=((B/A)^2+(C/A)^2*(sin_theta^2*cos_theta^2+sin_theta^4*cos_phi^2*sin_phi^2))^0.5;
        ei=(h*h*abs(A)/(2*m0))*sk*(1/q)*(1-g1);
        ef=ei;
        particle(1,7)=ei;
    elseif valley==4
        cos_theta=kz/ki;
        sin_theta=sqrt(1-cos_theta^2);
        cos_phi=kx/(ki*sin_theta);
        sin_phi=ky/(ki*sin_theta);
        g1=((B/A)^2+(C/A)^2*(sin_theta^2*cos_theta^2+sin_theta^4*cos_phi^2*sin_phi^2))^0.5;
        ei=(h*h*abs(A)/(2*m0))*sk*(1/q)*(1+g1);
        ef=ei;
        particle(1,7)=ei;
    end
    
    ie=abs(floor(ei/de)+1);
    if ie > iemax
        ie=iemax;
    end
    
    r1=rand();
    %------------------------------
    %Gamma Valley
    if valley==1
        %Acoustic Scattering
        if r1 < scatGaAs(1,ie,valley)
            kf=ki;
            phi=2*pi*rand();
            cos_t=1-2*rand();
            sin_t=sqrt(1-cos_t*cos_t);
            kxf=kf*sin_t*cos(phi);
            kyf=kf*sin_t*sin(phi);
            kzf=kf*cos_t;
%             valleyf=1;
        %POP Absorption
        elseif r1 < scatGaAs(2,ie,valley)
            ef=ei+hw0; %---------------------hwo or hwoq????????????
            ff=2*sqrt(ei*ef)/((sqrt(ei)-sqrt(ef))^2);
            
            if ff > 0
                kf=(sqrt(2*eM(valley))/h)*sqrt(ef*(1+alpha(valley)*ef))*sqrt(q);%Need to check if q is correct
                phi=2*pi*rand();
                cos_t=(1+ff-(1+2*ff)^rand())/ff;
                sin_t=sqrt(1-cos_t*cos_t);
            
                sin_a=sqrt(skx+sky)/ki;
                cos_a=kz/ki;
                sin_b=kx/sqrt(skx+sky);
                cos_b=ky/sqrt(skx+sky);
            
                kxr=kf*sin_t*cos(phi);
                kyr=kf*sin_t*sin(phi);
                kzr=kf*cos_t;
            
                kxf=cos_b*kxr+cos_a*sin_b*kyr+sin_a*sin_b*kzr;
                kyf=-sin_b*kxr+cos_a*cos_b*kyr+sin_a*cos_b*kzr;
                kzf=-sin_a*kyr+cos_a*kzr;
            end
%             valleyf=1;
            
        %POP Emission
        elseif r1 < scatGaAs(3,ie,valley)
            ef=ei-hw0; %---------------------hwo or hwoq????????????
            if ef > 0                
                ff=2*sqrt(ei*ef)/((sqrt(ei)-sqrt(ef))^2);
                
                if ff > 0
                    kf=(sqrt(2*eM(valley))/h)*sqrt(ef*(1+alpha(valley)*ef))*sqrt(q);%Need to check if q is correct
                    phi=2*pi*rand();
                    cos_t=(1+ff-(1+2*ff)^rand())/ff;
                    sin_t=sqrt(1-cos_t*cos_t);
            
                    sin_a=sqrt(skx+sky)/ki;
                    cos_a=kz/ki;
                    sin_b=kx/sqrt(skx+sky);
                    cos_b=ky/sqrt(skx+sky);
            
                    kxr=kf*sin_t*cos(phi);
                    kyr=kf*sin_t*sin(phi);
                    kzr=kf*cos_t;
            
                    kxf=cos_b*kxr+cos_a*sin_b*kyr+sin_a*sin_b*kzr;
                    kyf=-sin_b*kxr+cos_a*cos_b*kyr+sin_a*cos_b*kzr;
                    kzf=-sin_a*kyr+cos_a*kzr;
                end
            end
%             valleyf=1;
            
        %NPOP Absorption
        elseif r1 < scatGaAs(4,ie,valley)
            ef=ei+hwij-(Egl-Egg); %---------------------hwij or hwijq????????????
            valleyf=2;
            kf=(sqrt(2*eM(valleyf))/h)*sqrt(ef*(1+alpha(valleyf)*ef))*sqrt(q);%Need to check if q is correct
            phi=2*pi*rand();
            cos_t=1-2*rand();
            sin_t=sqrt(1-cos_t*cos_t);
            kxf=kf*sin_t*cos(phi);
            kyf=kf*sin_t*sin(phi);
            kzf=kf*cos_t;
        %NPOP Emission
        elseif r1 < scatGaAs(5,ie,valley)
            ef=ei-hwij-(Egl-Egg); %---------------------hwij or hwijq????????????
            valleyf=valley;
            if ef > 0
                valleyf=2;
                kf=(sqrt(2*eM(valleyf))/h)*sqrt(ef*(1+alpha(valleyf)*ef))*sqrt(q);%Need to check if q is correct
                phi=2*pi*rand();
                cos_t=1-2*rand();
                sin_t=sqrt(1-cos_t*cos_t);
                kxf=kf*sin_t*cos(phi);
                kyf=kf*sin_t*sin(phi);
                kzf=kf*cos_t;
            end
        %Impurity 
        elseif r1 < scatGaAs(6,ie,valley)
            kf=ki;
            r=rand();
            phi=2*pi*rand();
            cos_t=1-(2*r)/(1+(1-r)*(2*ki/qD)^2);
            sin_t=sqrt(1-cos_t*cos_t);

            sin_a=sqrt(skx+sky)/ki;
            cos_a=kz/ki;
            sin_b=kx/sqrt(skx+sky);
            cos_b=ky/sqrt(skx+sky);

            kxr=kf*sin_t*cos(phi);
            kyr=kf*sin_t*sin(phi);
            kzr=kf*cos_t;

            kxf=cos_b*kxr+cos_a*sin_b*kyr+sin_a*sin_b*kzr;
            kyf=-sin_b*kxr+cos_a*cos_b*kyr+sin_a*cos_b*kzr;
            kzf=-sin_a*kyr+cos_a*kzr;

            valleyf=1;
        end
        
    %L Valley--------------------------------------------------------------
    elseif valley==2
        %Acoustic
        if r1 < scatGaAs(1,ie,valley)
            kf=ki;
            phi=2*pi*rand();
            cos_t=1-2*rand();
            sin_t=sqrt(1-cos_t*cos_t);
            kxf=kf*sin_t*cos(phi);
            kyf=kf*sin_t*sin(phi);
            kzf=kf*cos_t;
            valleyf=2;
        %NPOP Absorption (Inter)
        elseif r1 < scatGaAs(2,ie,valley)
            ef=ei+hwij+(Egl-Egg); %---------------------hwij or hwijq????????????
            valleyf=1;
            kf=(sqrt(2*eM(valleyf))/h)*sqrt(ef*(1+alpha(valleyf)*ef))*sqrt(q);%Need to check if q is correct
            phi=2*pi*rand();
            cos_t=1-2*rand();
            sin_t=sqrt(1-cos_t*cos_t);
            kxf=kf*sin_t*cos(phi);
            kyf=kf*sin_t*sin(phi);
            kzf=kf*cos_t;
        %NPOP Emission (Inter)
        elseif r1 < scatGaAs(3,ie,valley)
            ef=ei-hwij+(Egl-Egg); %---------------------hwij or hwijq????????????
            valleyf=valley;
            if ef > 0
                valleyf=1;
                kf=(sqrt(2*eM(valleyf))/h)*sqrt(ef*(1+alpha(valleyf)*ef))*sqrt(q);%Need to check if q is correct
                phi=2*pi*rand();
                cos_t=1-2*rand();
                sin_t=sqrt(1-cos_t*cos_t);
                kxf=kf*sin_t*cos(phi);
                kyf=kf*sin_t*sin(phi);
                kzf=kf*cos_t;
            end
        %POP Absorption
        elseif r1 < scatGaAs(4,ie,valley)
            ef=ei+hw0; %---------------------hwo or hwoq????????????            
            ff=2*sqrt(ei*ef)/((sqrt(ei)-sqrt(ef))^2);
            
            if ff > 0
                kf=(sqrt(2*eM(valley))/h)*sqrt(ef*(1+alpha(valley)*ef))*sqrt(q);%Need to check if q is correct
                phi=2*pi*rand();
                cos_t=(1+ff-(1+2*ff)^rand())/ff;
                sin_t=sqrt(1-cos_t*cos_t);
            
                sin_a=sqrt(skx+sky)/ki;
                cos_a=kz/ki;
                sin_b=kx/sqrt(skx+sky);
                cos_b=ky/sqrt(skx+sky);
            
                kxr=kf*sin_t*cos(phi);
                kyr=kf*sin_t*sin(phi);
                kzr=kf*cos_t;
            
                kxf=cos_b*kxr+cos_a*sin_b*kyr+sin_a*sin_b*kzr;
                kyf=-sin_b*kxr+cos_a*cos_b*kyr+sin_a*cos_b*kzr;
                kzf=-sin_a*kyr+cos_a*kzr;
            end
            valleyf=2;
            
        %POP Emission
        elseif r1 < scatGaAs(5,ie,valley)
            ef=ei-hw0; %---------------------hwo or hwoq????????????
            if ef > 0                
                ff=2*sqrt(ei*ef)/((sqrt(ei)-sqrt(ef))^2);
                
                if ff > 0
                    kf=(sqrt(2*eM(valley))/h)*sqrt(ef*(1+alpha(valley)*ef))*sqrt(q);%Need to check if q is correct
                    phi=2*pi*rand();
                    cos_t=(1+ff-(1+2*ff)^rand())/ff;
                    sin_t=sqrt(1-cos_t*cos_t);
            
                    sin_a=sqrt(skx+sky)/ki;
                    cos_a=kz/ki;
                    sin_b=kx/sqrt(skx+sky);
                    cos_b=ky/sqrt(skx+sky);
            
                    kxr=kf*sin_t*cos(phi);
                    kyr=kf*sin_t*sin(phi);
                    kzr=kf*cos_t;
            
                    kxf=cos_b*kxr+cos_a*sin_b*kyr+sin_a*sin_b*kzr;
                    kyf=-sin_b*kxr+cos_a*cos_b*kyr+sin_a*cos_b*kzr;
                    kzf=-sin_a*kyr+cos_a*kzr;
                end
            end
            valleyf=2;
        %NPOP Absorption (Intra)
        elseif r1 < scatGaAs(6,ie,valley)
            ef=ei+hwe; %---------------------hwe or hweq????????????
            valleyf=2;
            kf=(sqrt(2*eM(valleyf))/h)*sqrt(ef*(1+alpha(valleyf)*ef))*sqrt(q);%Need to check if q is correct
            phi=2*pi*rand();
            cos_t=1-2*rand();
            sin_t=sqrt(1-cos_t*cos_t);
            kxf=kf*sin_t*cos(phi);
            kyf=kf*sin_t*sin(phi);
            kzf=kf*cos_t;
        %NPOP Emission (Intra)
        elseif r1 < scatGaAs(7,ie,valley)
            ef=ei-hwe; %---------------------hwe or hweq????????????
            valleyf=2;
            kf=(sqrt(2*eM(valleyf))/h)*sqrt(ef*(1+alpha(valleyf)*ef))*sqrt(q);%Need to check if q is correct
            phi=2*pi*rand();
            cos_t=1-2*rand();
            sin_t=sqrt(1-cos_t*cos_t);
            kxf=kf*sin_t*cos(phi);
            kyf=kf*sin_t*sin(phi);
            kzf=kf*cos_t;
        %Impurity
        elseif r1 < scatGaAs(8,ie,valley)
%             if x < Ll || x > (Ll+Lg)
                kf=ki;
                r=rand();
                phi=2*pi*rand();
                cos_t=1-(2*r)/(1+(1-r)*(2*ki/qD)^2);
                sin_t=sqrt(1-cos_t*cos_t);
            
                sin_a=sqrt(skx+sky)/ki;
                cos_a=kz/ki;
                sin_b=kx/sqrt(skx+sky);
                cos_b=ky/sqrt(skx+sky);
            
                kxr=kf*sin_t*cos(phi);
                kyr=kf*sin_t*sin(phi);
                kzr=kf*cos_t;
            
                kxf=cos_b*kxr+cos_a*sin_b*kyr+sin_a*sin_b*kzr;
                kyf=-sin_b*kxr+cos_a*cos_b*kyr+sin_a*cos_b*kzr;
                kzf=-sin_a*kyr+cos_a*kzr;
            
                valleyf=2;
        end
    %-----------------------------
    elseif valley==3 %heavy hole
        valley=valley-2;
        %impurity scattering, interband,don't have reference now
        if r1 < scatGaAs_hole(1,ie,valley)
            ef=ei;            
            phi=2*pi*rand();
            cos_t=1-2*rand();
            sin_t=sqrt(1-cos_t*cos_t);
            g2=((B/A)^2+(C/A)^2*(sin_t^2*cos_t^2+sin_t^4*cos(phi)^2*sin(phi)^2))^0.5;
            kf=sqrt(2*m0*q*ef/h/h/abs(A)/(1+g2));
            kxf=kf*sin_t*cos(phi);
            kyf=kf*sin_t*sin(phi);
            kzf=kf*cos_t;
            valleyf=4;
        elseif r1 < scatGaAs_hole(2,ie,valley)%impurity scattering, intraband, don't know 
            kf=ki;
            phi=2*pi*rand();
            cos_t=1-2*rand();
            sin_t=sqrt(1-cos_t*cos_t);
            kxf=kf*sin_t*cos(phi);
            kyf=kf*sin_t*sin(phi);
            kzf=kf*cos_t;
            valleyf=3;
        elseif r1 < scatGaAs_hole(3,ie,valley)%optical, h-h, absorbtion  
            r=rand();    
            ef=ei+hw0;            
            phi=2*pi*rand();
            cos_t=1-2*rand();
            sin_t=sqrt(1-cos_t*cos_t);
            g2=((B/A)^2+(C/A)^2*(sin_t^2*cos_t^2+sin_t^4*cos(phi)^2*sin(phi)^2))^0.5;
            kf=sqrt(2*m0*q*ef/h/h/abs(A)/(1-g2));
            kxf=kf*sin_t*cos(phi);
            kyf=kf*sin_t*sin(phi);
            kzf=kf*cos_t;          
            cos_alpha=(kx*kxf+ky*kyf+kz*kzf)/ki/kf;
            G=0.25*(1+3*cos_alpha^2);
            fa=((1-g2)/(1-g1)+ef/ei-2*cos_theta*sqrt((1-g2)/(1-g1)*ef/ei))*G;            
            fb=sqrt(1-g2)*((1-g2)/(1-g1)+ef/ei-2*cos_theta*sqrt((1-g2)/(1-g1)*ef/ei)+a*b*b*(1-g2))^2;
            f=fa/fb;
            fmax=(1+ef/ei-2*sqrt(ef/ei)*cos_theta)/(1-g111)^0.5/(1+ef/ei-2*sqrt(ef/ei)*cos_theta+a*b*b*(1-g111))^2;
            valleyf=3;
            i=1;
            while r>f/fmax
                phi=2*pi*rand();
                cos_t=1-2*rand();
                sin_t=sqrt(1-cos_t*cos_t);
                g2=((B/A)^2+(C/A)^2*(sin_t^2*cos_t^2+sin_t^4*cos(phi)^2*sin(phi)^2))^0.5;
                kf=sqrt(2*m0*q*ef/h/h/abs(A)/(1-g2));
                kxf=kf*sin_t*cos(phi);
                kyf=kf*sin_t*sin(phi);
                kzf=kf*cos_t;        
                cos_alpha=(kx*kxf+ky*kyf+kz*kzf)/ki/kf;
                G=0.25*(1+3*cos_alpha^2);
                fa=((1-g2)/(1-g1)+ef/ei-2*cos_theta*sqrt((1-g2)/(1-g1)*ef/ei))*G;            
                fb=sqrt(1-g2)*((1-g2)/(1-g1)+ef/ei-2*cos_theta*sqrt((1-g2)/(1-g1)*ef/ei)+a*b*b*(1-g2))^2;
                f=fa/fb;
                i=i+1;
                if i>2
                    kxf=kx;
                    kyf=ky;
                    kzf=kz;
                    ef=ei;
                    break
                end
            end            
        elseif r1 < scatGaAs_hole(4,ie,valley)%optical, h-l, absorbtion  
            r=rand();
            ef=ei+hw0;            
            phi=2*pi*rand();
            cos_t=1-2*rand();
            sin_t=sqrt(1-cos_t*cos_t);
            g2=((B/A)^2+(C/A)^2*(sin_t^2*cos_t^2+sin_t^4*cos(phi)^2*sin(phi)^2))^0.5;
            kf=sqrt(2*m0*q*ef/h/h/abs(A)/(1+g2));
            kxf=kf*sin_t*cos(phi);
            kyf=kf*sin_t*sin(phi);
            kzf=kf*cos_t;
            cos_alpha=(kx*kxf+ky*kyf+kz*kzf)/ki/kf;
            G=0.75*(1-cos_alpha^2);
            fa=((1+g2)/(1-g1)+ef/ei-2*cos_theta*sqrt((1+g2)/(1-g1)*ef/ei))*G;
            fb=sqrt(1+g2)*((1+g2)/(1-g1)+ef/ei-2*cos_theta*sqrt((1+g2)/(1-g1)*ef/ei)+a*b*b*(1+g2))^2;
            f=fa/fb;
            fmaxa=0.75*((1+g111)/(1-g111)+(ef/ei)-2*cos_theta*sqrt((1+g111)/(1-g111)*ef/ei));
            fmaxb=sqrt(1+g100)*((1+g100)/(1-g100)+(ef/ei)-2*cos_theta*sqrt((1+g100)/(1-g100)*ef/ei)+a*b*b*(1+g100))^2;
            fmax=fmaxa/fmaxb;
            valleyf=4;
            i=1;
            while r>f/fmax
                phi=2*pi*rand();
                cos_t=1-2*rand();
                sin_t=sqrt(1-cos_t*cos_t);
                g2=((B/A)^2+(C/A)^2*(sin_t^2*cos_t^2+sin_t^4*cos(phi)^2*sin(phi)^2))^0.5;
                kf=sqrt(2*m0*q*ef/h/h/abs(A)/(1+g2));
                kxf=kf*sin_t*cos(phi);
                kyf=kf*sin_t*sin(phi);
                kzf=kf*cos_t;
                cos_alpha=(kx*kxf+ky*kyf+kz*kzf)/ki/kf;
                G=0.75*(1-cos_alpha^2);
                fa=((1+g2)/(1-g1)+ef/ei-2*cos_theta*sqrt((1+g2)/(1-g1)*ef/ei))*G;
                fb=sqrt(1+g2)*((1+g2)/(1-g1)+ef/ei-2*cos_theta*sqrt((1+g2)/(1-g1)*ef/ei)+a*b*b*(1+g2))^2;
                f=fa/fb;
                i=i+1;
                if i>2
                    kxf=kx;
                    kyf=ky;
                    kzf=kz;
                    ef=ei;
                    valleyf=3;
                    break
                end
            end            
        elseif r1 < scatGaAs_hole(5,ie,valley)%optical, h-h, emission
            r=rand();         
            ef=ei-hw0; 
            if ef<0
                valleyf=3;%self scattering
                fprintf('5h is wrong,ef=%f\n',ef);%914
            else
                phi=2*pi*rand();
                cos_t=1-2*rand();
                sin_t=sqrt(1-cos_t*cos_t);
                g2=((B/A)^2+(C/A)^2*(sin_t^2*cos_t^2+sin_t^4*cos(phi)^2*sin(phi)^2))^0.5;
                kf=sqrt(2*m0*q*ef/h/h/abs(A)/(1-g2));
                kxf=kf*sin_t*cos(phi);
                kyf=kf*sin_t*sin(phi);
                kzf=kf*cos_t;
                g2=((B/A)^2+(C/A)^2*(sin_t^2*cos_t^2+sin_t^4*cos(phi)^2*sin(phi)^2))^0.5;
                cos_alpha=(kx*kxf+ky*kyf+kz*kzf)/ki/kf;
                G=0.25*(1+3*cos_alpha^2);
                fa=((1-g2)/(1-g1)+ef/ei-2*cos_theta*sqrt((1-g2)/(1-g1)*ef/ei))*G;            
                fb=sqrt(1-g2)*((1-g2)/(1-g1)+ef/ei-2*cos_theta*sqrt((1-g2)/(1-g1)*ef/ei)+a*b*b*(1-g2))^2;
                f=fa/fb;
                valleyf=3;
                fmax=(1+ef/ei-2*sqrt(ef/ei)*cos_theta)/(1-g111)^0.5/(1+ef/ei-2*sqrt(ef/ei)*cos_theta+a*b*b*(1-g111))^2;
                i=1;
                while r>f/fmax
                    phi=2*pi*rand();
                    cos_t=1-2*rand();
                    sin_t=sqrt(1-cos_t*cos_t);
                    g2=((B/A)^2+(C/A)^2*(sin_t^2*cos_t^2+sin_t^4*cos(phi)^2*sin(phi)^2))^0.5;
                    kf=sqrt(2*m0*q*ef/h/h/abs(A)/(1-g2));
                    kxf=kf*sin_t*cos(phi);
                    kyf=kf*sin_t*sin(phi);
                    kzf=kf*cos_t;
                    g2=((B/A)^2+(C/A)^2*(sin_t^2*cos_t^2+sin_t^4*cos(phi)^2*sin(phi)^2))^0.5;
                    cos_alpha=(kx*kxf+ky*kyf+kz*kzf)/ki/kf;
                    G=0.25*(1+3*cos_alpha^2);
                    fa=((1-g2)/(1-g1)+ef/ei-2*cos_theta*sqrt((1-g2)/(1-g1)*ef/ei))*G;            
                    fb=sqrt(1-g2)*((1-g2)/(1-g1)+ef/ei-2*cos_theta*sqrt((1-g2)/(1-g1)*ef/ei)+a*b*b*(1-g2))^2;
                    f=fa/fb;
                    i=i+1;
                    if i>2
                        kxf=kx;
                        kyf=ky;
                        kzf=kz;
                        ef=ei;
                        break
                    end
                end
            end            
        elseif r1 < scatGaAs_hole(6,ie,valley)%optical, h-l, emission   
            r=rand();
            ef=ei-hw0;
            if ef<0
                fprintf('6h is wrong,ef=%f\n',ef);%914
                ef=ef+de;
                valleyf=4;%self scattering
            else
                phi=2*pi*rand();
                cos_t=1-2*rand();
                sin_t=sqrt(1-cos_t*cos_t);
                g2=((B/A)^2+(C/A)^2*(sin_t^2*cos_t^2+sin_t^4*cos(phi)^2*sin(phi)^2))^0.5;
                kf=sqrt(2*m0*q*ef/h/h/abs(A)/(1+g2));
                kxf=kf*sin_t*cos(phi);
                kyf=kf*sin_t*sin(phi);
                kzf=kf*cos_t;
                g2=((B/A)^2+(C/A)^2*(sin_t^2*cos_t^2+sin_t^4*cos(phi)^2*sin(phi)^2))^0.5;
                cos_alpha=(kx*kxf+ky*kyf+kz*kzf)/ki/kf;
                G=0.75*(1-cos_alpha^2);
                fa=((1+g2)/(1-g1)+ef/ei-2*cos_theta*sqrt((1+g2)/(1-g1)*ef/ei))*G;
                fb=sqrt(1+g2)*((1+g2)/(1-g1)+ef/ei-2*cos_theta*sqrt((1+g2)/(1-g1)*ef/ei)+a*b*b*(1+g2))^2;
                f=fa/fb;
                g111=sqrt((B/A)^2+(C/A)^2/3);
                g100=B/A;
                fmaxa=0.75*((1+g111)/(1-g111)+(ef/ei)-2*cos_theta*sqrt((1+g111)/(1-g111)*ef/ei));
                fmaxb=sqrt(1+g100)*((1+g100)/(1-g100)+(ef/ei)-2*cos_theta*sqrt((1+g100)/(1-g100)*ef/ei)+a*b*b*(1+g100))^2;
                fmax=fmaxa/fmaxb;
                valleyf=4;
                i=1;
                while r>f/fmax
                    phi=2*pi*rand();
                    cos_t=1-2*rand();
                    sin_t=sqrt(1-cos_t*cos_t);
                    g2=((B/A)^2+(C/A)^2*(sin_t^2*cos_t^2+sin_t^4*cos(phi)^2*sin(phi)^2))^0.5;
                    kf=sqrt(2*m0*q*ef/h/h/abs(A)/(1+g2));
                    kxf=kf*sin_t*cos(phi);
                    kyf=kf*sin_t*sin(phi);
                    kzf=kf*cos_t;
                    g2=((B/A)^2+(C/A)^2*(sin_t^2*cos_t^2+sin_t^4*cos(phi)^2*sin(phi)^2))^0.5;
                    cos_alpha=(kx*kxf+ky*kyf+kz*kzf)/ki/kf;
                    G=0.75*(1-cos_alpha^2);
                    fa=((1+g2)/(1-g1)+ef/ei-2*cos_theta*sqrt((1+g2)/(1-g1)*ef/ei))*G;
                    fb=sqrt(1+g2)*((1+g2)/(1-g1)+ef/ei-2*cos_theta*sqrt((1+g2)/(1-g1)*ef/ei)+a*b*b*(1+g2))^2;
                    f=fa/fb;
                    i=i+1;
                    if i>2
                        kxf=kx;
                        kyf=ky;
                        kzf=kz;
                        valleyf=3;
                        ef=ei;
                        break
                    end
                end
            end
        elseif r1 < scatGaAs_hole(7,ie,valley)%non-polar, to h,absorb
            r=rand();
            flag=0;
            phi=2*pi*rand();
            cos_t=1-2*rand();
            sin_t=sqrt(1-cos_t*cos_t);
            g=((B/A)^2+(C/A)^2*(sin_t*sin_t*cos_t*cos_t+sin_t^4*cos(phi)^2*sin(phi)^2))^0.5;
            fmax=1/(a*(1-sqrt((B/A)^2+(C/A)^2/3)))^1.5;
            f=1/(a*(1-g))^1.5;
            i=1;
            while r>f/fmax
                phi=2*pi*rand();
                cos_t=1-2*rand();
                sin_t=sqrt(1-cos_t*cos_t);
                g=((B/A)^2+(C/A)^2*(sin_t*sin_t*cos_t*cos_t+sin_t^4*cos(phi)^2*sin(phi)^2))^0.5;
                f=1/(a*(1-g))^1.5;
                i=i+1;
                if i>2
                    flag=1;
                    break;
                end
            end
            if flag==0
                ef=ei+hw0;
%                 g=((B/A)^2+(C/A)^2*(sin_t^2*cos_t^2+sin_t^4*cos(phi)^2*sin(phi)^2))^0.5;
                kf=sqrt(ef*2*m0*q/abs(A)/(1-g))/h;
                kxf=kf*sin_t*cos(phi);
                kyf=kf*sin_t*sin(phi);
                kzf=kf*cos_t;
                valleyf=3;
            else
                valleyf=3;
            end
        elseif r1 < scatGaAs_hole(8,ie,valley) %non-polar, to h,emisson
            ef=ei-hw0;
            if ef<0
                valleyf=3;
            else
                r=rand();
                flag=0;
                phi=2*pi*rand();
                cos_t=1-2*rand();
                sin_t=sqrt(1-cos_t*cos_t);
                g=((B/A)^2+(C/A)^2*(sin_t*sin_t*cos_t*cos_t+sin_t^4*cos(phi)^2*sin(phi)^2))^0.5;
                fmax=1/(a*(1-sqrt((B/A)^2+(C/A)^2/3)))^1.5;
                f=1/(a*(1-g))^1.5;
                i=1;
                while r>f/fmax
                    phi=2*pi*rand();
                    cos_t=1-2*rand();
                    sin_t=sqrt(1-cos_t*cos_t);
                    g=((B/A)^2+(C/A)^2*(sin_t*sin_t*cos_t*cos_t+sin_t^4*cos(phi)^2*sin(phi)^2))^0.5;
                    f=1/(a*(1-g))^1.5;
                    i=i+1;
                    if i>2
                        flag=1;
                        break;
                    end
                end
                if flag==0
%                     g=((B/A)^2+(C/A)^2*(sin_t^2*cos_t^2+sin_t^4*cos(phi)^2*sin(phi)^2))^0.5;
                    kf=sqrt(ef*2*m0*q/abs(A)/(1-g))/h;
                    kxf=kf*sin_t*cos(phi);
                    kyf=kf*sin_t*sin(phi);
                    kzf=kf*cos_t;
                    valleyf=3;
                else
                    valleyf=3;
                end
            end          
        elseif r1 < scatGaAs_hole(9,ie,valley)%non-polar, to l, aorbsion
            r=rand();
            flag=0;
            phi=2*pi*rand();
            cos_t=1-2*rand();
            sin_t=sqrt(1-cos_t*cos_t);
            g=((B/A)^2+(C/A)^2*(sin_t*sin_t*cos_t*cos_t+sin_t^4*cos(phi)^2*sin(phi)^2))^0.5;
            fmax=1/(a*(1+B/A))^1.5;
            f=1/(a*(1+g))^1.5;
            i=1;
            while r>f/fmax
                phi=2*pi*rand();
                cos_t=1-2*rand();
                sin_t=sqrt(1-cos_t*cos_t);
                g=((B/A)^2+(C/A)^2*(sin_t*sin_t*cos_t*cos_t+sin_t^4*cos(phi)^2*sin(phi)^2))^0.5;
                f=1/(a*(1-g))^1.5;
                i=i+1;
                if i>2
                    flag=1;
                    break;
                end
            end
            if flag==0
                ef=ei+hw0;
%                 g=((B/A)^2+(C/A)^2*(sin_t^2*cos_t^2+sin_t^4*cos(phi)^2*sin(phi)^2))^0.5;
                kf=sqrt(ef*2*m0*q/abs(A)/(1+g))/h;
                kxf=kf*sin_t*cos(phi);
                kyf=kf*sin_t*sin(phi);
                kzf=kf*cos_t;
                valleyf=4;
            else
                valleyf=3;
            end
        elseif r1 < scatGaAs_hole(10,ie,valley)%non-polar, to l, emission
            ef=ei-hw0;
            if ef<0
                 valleyf=3;
            else
                r=rand();
                flag=0;
                phi=2*pi*rand();
                cos_t=1-2*rand();
                sin_t=sqrt(1-cos_t*cos_t);
                g=((B/A)^2+(C/A)^2*(sin_t*sin_t*cos_t*cos_t+sin_t^4*cos(phi)^2*sin(phi)^2))^0.5;
                fmax=1/(a*(1+B/A))^1.5;
                f=1/(a*(1+g))^1.5;
                i=1;
                while r>f/fmax
                    phi=2*pi*rand();
                    cos_t=1-2*rand();
                    sin_t=sqrt(1-cos_t*cos_t);
                    g=((B/A)^2+(C/A)^2*(sin_t*sin_t*cos_t*cos_t+sin_t^4*cos(phi)^2*sin(phi)^2))^0.5;
                    f=1/(a*(1-g))^1.5;
                    i=i+1;
                    if i>2
                        flag=1;
                        break;
                    end
                end
                if flag==0
%                     g=((B/A)^2+(C/A)^2*(sin_t^2*cos_t^2+sin_t^4*cos(phi)^2*sin(phi)^2))^0.5;
                    kf=sqrt(ef*2*m0*q/abs(A)/(1+g))/h;
                    kxf=kf*sin_t*cos(phi);
                    kyf=kf*sin_t*sin(phi);
                    kzf=kf*cos_t;
                    valleyf=4;
                else
                    valleyf=3;
                end
            end          
        elseif r1 < scatGaAs_hole(11,ie,valley)%acoustic, intra, h-h
            r=rand();
            flag=0;
            phi=2*pi*rand();
            cos_t=1-2*rand();
            sin_t=sqrt(1-cos_t*cos_t);
            g=((B/A)^2+(C/A)^2*(sin_t*sin_t*cos_t*cos_t+sin_t^4*cos(phi)^2*sin(phi)^2))^0.5;
            fmax=1/(a*(1-sqrt((B/A)^2+(C/A)^2/3)))^1.5;
            f=1/(a*(1-g))^1.5;
            i=1;
            while r>f/fmax
                phi=2*pi*rand();
                cos_t=1-2*rand();
                sin_t=sqrt(1-cos_t*cos_t);
                g=((B/A)^2+(C/A)^2*(sin_t*sin_t*cos_t*cos_t+sin_t^4*cos(phi)^2*sin(phi)^2))^0.5;
                f=1/(a*(1-g))^1.5;
                i=i+1;
                if i>2
                    flag=1;
                    break;
                end
            end
            if flag==0
                ef=ei;
%                 g=((B/A)^2+(C/A)^2*(sin_t^2*cos_t^2+sin_t^4*cos(phi)^2*sin(phi)^2))^0.5;
                kf=sqrt(ef*2*m0*q/abs(A)/(1-g))/h;
                kxf=kf*sin_t*cos(phi);
                kyf=kf*sin_t*sin(phi);
                kzf=kf*cos_t;
                valleyf=3;
            else
                valleyf=3;
            end
        elseif r1 < scatGaAs_hole(12,ie,valley)%acoustic, inter,h-l
            r=rand();
            flag=0;
            phi=2*pi*rand();
            cos_t=1-2*rand();
            sin_t=sqrt(1-cos_t*cos_t);
            g=((B/A)^2+(C/A)^2*(sin_t*sin_t*cos_t*cos_t+sin_t^4*cos(phi)^2*sin(phi)^2))^0.5;
            fmax=1/(a*(1+B/A))^1.5;
            f=1/(a*(1+g))^1.5;
            i=1;
            while r>f/fmax
                phi=2*pi*rand();
                cos_t=1-2*rand();
                sin_t=sqrt(1-cos_t*cos_t);
                g=((B/A)^2+(C/A)^2*(sin_t*sin_t*cos_t*cos_t+sin_t^4*cos(phi)^2*sin(phi)^2))^0.5;
                f=1/(a*(1-g))^1.5;
                i=i+1;
                if i>2
                    flag=1;
                    break;
                end
            end
            if flag==0
                ef=ei;
%                 g=((B/A)^2+(C/A)^2*(sin_t^2*cos_t^2+sin_t^4*cos(phi)^2*sin(phi)^2))^0.5;
                kf=sqrt(ef*2*m0*q/abs(A)/(1+g))/h;
                kxf=kf*sin_t*cos(phi);
                kyf=kf*sin_t*sin(phi);
                kzf=kf*cos_t;
                valleyf=4;
            else
                valleyf=3;
            end
        end
            
        %-----------------------------------
    elseif valley==4 %light hole
        valley=valley-2;
        %impurity scattering, l-h
        if r1 < scatGaAs_hole(1,ie,valley)
            ef=ei;
            phi=2*pi*rand();
            cos_t=1-2*rand();
            sin_t=sqrt(1-cos_t*cos_t);
            g=((B/A)^2+(C/A)^2*(sin_t^2*cos_t^2+sin_t^4*cos(phi)^2*sin(phi)^2))^0.5;
            kf=sqrt(ef*2*m0*q/abs(A)/(1-g))/h;
            kxf=kf*sin_t*cos(phi);
            kyf=kf*sin_t*sin(phi);
            kzf=kf*cos_t;
            valleyf=3;
        elseif r1 < scatGaAs_hole(2,ie,valley)%impurity scattering, l-l
            ef=ei;            
            phi=2*pi*rand();
            cos_t=1-2*rand();
            sin_t=sqrt(1-cos_t*cos_t);
            g=((B/A)^2+(C/A)^2*(sin_t^2*cos_t^2+sin_t^4*cos(phi)^2*sin(phi)^2))^0.5;
            kf=sqrt(ef*2*m0*q/abs(A)/(1+g))/h;
            kxf=kf*sin_t*cos(phi);
            kyf=kf*sin_t*sin(phi);
            kzf=kf*cos_t;
            valleyf=4;
        elseif r1 < scatGaAs_hole(3,ie,valley)%optical, l-l, absorbtion  
            r=rand();
            ef=ei+hw0;           
            phi=2*pi*rand();
            cos_t=1-2*rand();
            sin_t=sqrt(1-cos_t*cos_t);
            g=((B/A)^2+(C/A)^2*(sin_t^2*cos_t^2+sin_t^4*cos(phi)^2*sin(phi)^2))^0.5;
            kf=sqrt(ef*2*m0*q/abs(A)/(1+g))/h;
            kxf=kf*sin_t*cos(phi);
            kyf=kf*sin_t*sin(phi);
            kzf=kf*cos_t;
            g2=((B/A)^2+(C/A)^2*(sin_t^2*cos_t^2+sin_t^4*cos(phi)^2*sin(phi)^2))^0.5;
            cos_alpha=(kx*kxf+ky*kyf+kz*kzf)/ki/kf;
            G=0.25*(1+3*cos_alpha^2);
            fa=((1+g2)/(1+g1)+ef/ei-2*cos_theta*sqrt((1+g2)/(1+g1)*ef/ei))*G;
            fb=sqrt(1+g2)*((1+g2)/(1+g1)+ef/ei-2*cos_theta*sqrt((1+g2)/(1+g1)*ef/ei)+a*b*b*(1+g2))^2;
            f=fa/fb;
            fmax=(1+ef/ei-2*sqrt(ef/ei)*cos_theta)/(1+g100)^0.5/(1+ef/ei-2*sqrt(ef/ei)*cos_theta+a*b*b*(1+g100))^2;
            valleyf=4;
            i=1;
            while r>f/fmax
                phi=2*pi*rand();
                cos_t=1-2*rand();
                sin_t=sqrt(1-cos_t*cos_t);
                g=((B/A)^2+(C/A)^2*(sin_t^2*cos_t^2+sin_t^4*cos(phi)^2*sin(phi)^2))^0.5;
                kf=sqrt(ef*2*m0*q/abs(A)/(1+g))/h;
                kxf=kf*sin_t*cos(phi);
                kyf=kf*sin_t*sin(phi);
                kzf=kf*cos_t;
                g2=((B/A)^2+(C/A)^2*(sin_t^2*cos_t^2+sin_t^4*cos(phi)^2*sin(phi)^2))^0.5;
                cos_alpha=(kx*kxf+ky*kyf+kz*kzf)/ki/kf;
                G=0.25*(1+3*cos_alpha^2);
                fa=((1+g2)/(1+g1)+ef/ei-2*cos_theta*sqrt((1+g2)/(1+g1)*ef/ei))*G;
                fb=sqrt(1+g2)*((1+g2)/(1+g1)+ef/ei-2*cos_theta*sqrt((1+g2)/(1+g1)*ef/ei)+a*b*b*(1+g2))^2;
                f=fa/fb;
                fmax=(1+ef/ei-2*sqrt(ef/ei)*cos_theta)/(1+g100)^0.5/(1+ef/ei-2*sqrt(ef/ei)*cos_theta+a*b*b*(1+g100))^2;
                valleyf=4;
                i=i+1;
                if i>2
                    kxf=kx;
                    kyf=ky;
                    kzf=kz;
                    break
                end
            end            
        elseif r1 < scatGaAs_hole(4,ie,valley)%optical, l-h, absorbtion  
            r=rand();
            ef=ei+hw0;            
            phi=2*pi*rand();
            cos_t=1-2*rand();
            sin_t=sqrt(1-cos_t*cos_t);
            g=((B/A)^2+(C/A)^2*(sin_t^2*cos_t^2+sin_t^4*cos(phi)^2*sin(phi)^2))^0.5;
            kf=sqrt(ef*2*m0*q/abs(A)/(1-g))/h;
            kxf=kf*sin_t*cos(phi);
            kyf=kf*sin_t*sin(phi);
            kzf=kf*cos_t;
            g2=((B/A)^2+(C/A)^2*(sin_t^2*cos_t^2+sin_t^4*cos(phi)^2*sin(phi)^2))^0.5;
            cos_alpha=(kx*kxf+ky*kyf+kz*kzf)/ki/kf;
            G=0.75*(1-cos_alpha^2);
            fa=((1-g2)/(1+g1)+ef/ei-2*cos_theta*sqrt((1-g2)/(1+g1)*ef/ei))*G;
            fb=sqrt(1-g2)*((1-g2)/(1+g1)+ef/ei-2*cos_theta*sqrt((1-g2)/(1+g1)*ef/ei)+a*b*b*(1-g2))^2;
            f=fa/fb;
            fmaxa=0.75*((1-g100)/(1+g100)+(ef/ei)-2*cos_theta*sqrt((1-g100)/(1+g100)*ef/ei));
            fmaxb=sqrt(1-g111)*((1-g111)/(1+g111)+(ef/ei)-2*cos_theta*sqrt((1-g111)/(1+g111)*ef/ei)+a*b*b*(1-g111))^2;
            fmax=fmaxa/fmaxb;
            i=1;
            valleyf=3;
            while r>f/fmax
                phi=2*pi*rand();
                cos_t=1-2*rand();
                sin_t=sqrt(1-cos_t*cos_t);
                g=((B/A)^2+(C/A)^2*(sin_t^2*cos_t^2+sin_t^4*cos(phi)^2*sin(phi)^2))^0.5;
                kf=sqrt(ef*2*m0*q/abs(A)/(1-g))/h;
                kxf=kf*sin_t*cos(phi);
                kyf=kf*sin_t*sin(phi);
                kzf=kf*cos_t;
                g2=((B/A)^2+(C/A)^2*(sin_t^2*cos_t^2+sin_t^4*cos(phi)^2*sin(phi)^2))^0.5;
                cos_alpha=(kx*kxf+ky*kyf+kz*kzf)/ki/kf;
                G=0.75*(1-cos_alpha^2);
                fa=((1-g2)/(1+g1)+ef/ei-2*cos_theta*sqrt((1-g2)/(1+g1)*ef/ei))*G;
                fb=sqrt(1-g2)*((1-g2)/(1+g1)+ef/ei-2*cos_theta*sqrt((1-g2)/(1+g1)*ef/ei)+a*b*b*(1-g2))^2;
                f=fa/fb;
                fmaxa=0.75*((1-g100)/(1+g100)+(ef/ei)-2*cos_theta*sqrt((1-g100)/(1+g100)*ef/ei));
                fmaxb=sqrt(1-g111)*((1-g111)/(1+g111)+(ef/ei)-2*cos_theta*sqrt((1-g111)/(1+g111)*ef/ei)+a*b*b*(1-g111))^2;
                fmax=fmaxa/fmaxb;
                valleyf=3;
                i=i+1;
                if i>2
                    kxf=kx;
                    kyf=ky;
                    kzf=kz;
                    valleyf=4;
                    break
                end
            end            
        elseif r1 < scatGaAs_hole(5,ie,valley)%optical, l-l, emission
            r=rand();
            ef=ei-hw0;
            if ef<0
                valleyf=4;
            else
                phi=2*pi*rand();
                cos_t=1-2*rand();
                sin_t=sqrt(1-cos_t*cos_t);
                g=((B/A)^2+(C/A)^2*(sin_t^2*cos_t^2+sin_t^4*cos(phi)^2*sin(phi)^2))^0.5;
                kf=sqrt(ef*2*m0*q/abs(A)/(1+g))/h;
                kxf=kf*sin_t*cos(phi);
                kyf=kf*sin_t*sin(phi);
                kzf=kf*cos_t;
                g2=((B/A)^2+(C/A)^2*(sin_t^2*cos_t^2+sin_t^4*cos(phi)^2*sin(phi)^2))^0.5;
                cos_alpha=(kx*kxf+ky*kyf+kz*kzf)/ki/kf;
                G=0.25*(1+3*cos_alpha^2);
                fa=((1+g2)/(1+g1)+ef/ei-2*cos_theta*sqrt((1+g2)/(1+g1)*ef/ei))*G;
                fb=sqrt(1+g2)*((1+g2)/(1+g1)+ef/ei-2*cos_theta*sqrt((1+g2)/(1+g1)*ef/ei)+a*b*b*(1+g2))^2;
                f=fa/fb;
                fmax=(1+ef/ei-2*sqrt(ef/ei)*cos_theta)/(1+g100)^0.5/(1+ef/ei-2*sqrt(ef/ei)*cos_theta+a*b*b*(1+g100))^2;
                i=1;
                valleyf=4;
                while r>f/fmax
                    phi=2*pi*rand();
                    cos_t=1-2*rand();
                    sin_t=sqrt(1-cos_t*cos_t);
                    g=((B/A)^2+(C/A)^2*(sin_t^2*cos_t^2+sin_t^4*cos(phi)^2*sin(phi)^2))^0.5;
                    kf=sqrt(ef*2*m0*q/abs(A)/(1+g))/h;
                    kxf=kf*sin_t*cos(phi);
                    kyf=kf*sin_t*sin(phi);
                    kzf=kf*cos_t;
                    g2=((B/A)^2+(C/A)^2*(sin_t^2*cos_t^2+sin_t^4*cos(phi)^2*sin(phi)^2))^0.5;
                    cos_alpha=(kx*kxf+ky*kyf+kz*kzf)/ki/kf;
                    G=0.25*(1+3*cos_alpha^2);
                    fa=((1+g2)/(1+g1)+ef/ei-2*cos_theta*sqrt((1+g2)/(1+g1)*ef/ei))*G;
                    fb=sqrt(1+g2)*((1+g2)/(1+g1)+ef/ei-2*cos_theta*sqrt((1+g2)/(1+g1)*ef/ei)+a*b*b*(1+g2))^2;
                    f=fa/fb;
                    fmax=(1+ef/ei-2*sqrt(ef/ei)*cos_theta)/(1+g100)^0.5/(1+ef/ei-2*sqrt(ef/ei)*cos_theta+a*b*b*(1+g100))^2;
                    valleyf=4;
                    i=i+1;
                    if i>2
                        kxf=kx;
                        kyf=ky;
                        kzf=kz;
                        valleyf=4;
                        break
                    end
                end
            end            
        elseif r1 < scatGaAs_hole(6,ie,valley)%optical, l-h, emission   
            r=rand();
            ef=ei-hw0;
            if ef<0
                valleyf=4;
            else                
                phi=2*pi*rand();
                cos_t=1-2*rand();
                sin_t=sqrt(1-cos_t*cos_t);
                g=((B/A)^2+(C/A)^2*(sin_t^2*cos_t^2+sin_t^4*cos(phi)^2*sin(phi)^2))^0.5;
                kf=sqrt(ef*2*m0*q/abs(A)/(1-g))/h;
                kxf=kf*sin_t*cos(phi);
                kyf=kf*sin_t*sin(phi);
                kzf=kf*cos_t;
                g2=((B/A)^2+(C/A)^2*(sin_t^2*cos_t^2+sin_t^4*cos(phi)^2*sin(phi)^2))^0.5;
                cos_alpha=(kx*kxf+ky*kyf+kz*kzf)/ki/kf;
                G=0.75*(1-cos_alpha^2);
                fa=((1-g2)/(1+g1)+ef/ei-2*cos_theta*sqrt((1-g2)/(1+g1)*ef/ei))*G;
                fb=sqrt(1-g2)*((1-g2)/(1+g1)+ef/ei-2*cos_theta*sqrt((1-g2)/(1+g1)*ef/ei)+a*b*b*(1-g2))^2;
                f=fa/fb;
                fmaxa=0.75*((1-g100)/(1+g100)+(ef/ei)-2*cos_theta*sqrt((1-g100)/(1+g100)*ef/ei));
                fmaxb=sqrt(1-g111)*((1-g111)/(1+g111)+(ef/ei)-2*cos_theta*sqrt((1-g111)/(1+g111)*ef/ei)+a*b*b*(1-g111))^2;
                fmax=fmaxa/fmaxb;
                i=1;
                valleyf=3;
                while r>f/fmax
                    phi=2*pi*rand();
                    cos_t=1-2*rand();
                    sin_t=sqrt(1-cos_t*cos_t);
                    g=((B/A)^2+(C/A)^2*(sin_t^2*cos_t^2+sin_t^4*cos(phi)^2*sin(phi)^2))^0.5;
                    kf=sqrt(ef*2*m0*q/abs(A)/(1-g))/h;
                    kxf=kf*sin_t*cos(phi);
                    kyf=kf*sin_t*sin(phi);
                    kzf=kf*cos_t;
                    g2=((B/A)^2+(C/A)^2*(sin_t^2*cos_t^2+sin_t^4*cos(phi)^2*sin(phi)^2))^0.5;
                    cos_alpha=(kx*kxf+ky*kyf+kz*kzf)/ki/kf;
                    G=0.75*(1-cos_alpha^2);
                    fa=((1-g2)/(1+g1)+ef/ei-2*cos_theta*sqrt((1-g2)/(1+g1)*ef/ei))*G;
                    fb=sqrt(1-g2)*((1-g2)/(1+g1)+ef/ei-2*cos_theta*sqrt((1-g2)/(1+g1)*ef/ei)+a*b*b*(1-g2))^2;
                    f=fa/fb;
                    fmaxa=0.75*((1-g100)/(1+g100)+(ef/ei)-2*cos_theta*sqrt((1-g100)/(1+g100)*ef/ei));
                    fmaxb=sqrt(1-g111)*((1-g111)/(1+g111)+(ef/ei)-2*cos_theta*sqrt((1-g111)/(1+g111)*ef/ei)+a*b*b*(1-g111))^2;
                    fmax=fmaxa/fmaxb;
                    valleyf=3;
                    i=i+1;
                    if i>2
                        kxf=kx;
                        kyf=ky;
                        kzf=kz;
                        valleyf=4;
                        break
                    end
                end
            end            
        elseif r1 < scatGaAs_hole(7,ie,valley)%non-polar, to h,absorb
            r=rand();
            flag=0;
            phi=2*pi*rand();
            cos_t=1-2*rand();
            sin_t=sqrt(1-cos_t*cos_t);
            g=((B/A)^2+(C/A)^2*(sin_t*sin_t*cos_t*cos_t+sin_t^4*cos(phi)^2*sin(phi)^2))^0.5;
            fmax=1/(a*(1+B/A))^1.5;
            f=1/(a*(1+g))^1.5;
            i=1;
            while r>f/fmax
                phi=2*pi*rand();
                cos_t=1-2*rand();
                sin_t=sqrt(1-cos_t*cos_t);
                g=((B/A)^2+(C/A)^2*(sin_t*sin_t*cos_t*cos_t+sin_t^4*cos(phi)^2*sin(phi)^2))^0.5;
                f=1/(a*(1+g))^1.5;
                i=i+1;
                if i>2
                    flag=1;
                    break;
                end
            end
            if flag==0
                ef=ei+hw0;
%                 g=((B/A)^2+(C/A)^2*(sin_t^2*cos_t^2+sin_t^4*cos(phi)^2*sin(phi)^2))^0.5;
                kf=sqrt(ef*2*m0*q/abs(A)/(1-g))/h;
                kxf=kf*sin_t*cos(phi);
                kyf=kf*sin_t*sin(phi);
                kzf=kf*cos_t;
                valleyf=3;     
            else
                valleyf=4;
            end                
        elseif r1 < scatGaAs_hole(8,ie,valley) %non-polar, to h,emisson
            r=rand();
            ef=ei-hw0;
            if ef<0
                valleyf=4;
            else
                flag=0;
                phi=2*pi*rand();
                cos_t=1-2*rand();
                sin_t=sqrt(1-cos_t*cos_t);
                g=((B/A)^2+(C/A)^2*(sin_t*sin_t*cos_t*cos_t+sin_t^4*cos(phi)^2*sin(phi)^2))^0.5;
                fmax=1/(a*(1+B/A))^1.5;
                f=1/(a*(1+g))^1.5;
                i=1;
                while r>f/fmax
                    phi=2*pi*rand();
                    cos_t=1-2*rand();
                    sin_t=sqrt(1-cos_t*cos_t);
                    g=((B/A)^2+(C/A)^2*(sin_t*sin_t*cos_t*cos_t+sin_t^4*cos(phi)^2*sin(phi)^2))^0.5;
                    f=1/(a*(1+g))^1.5;
                    i=i+1;
                    if i>2
                        flag=1;
                        break;
                    end
                end
                if flag==0
%                     g=((B/A)^2+(C/A)^2*(sin_t^2*cos_t^2+sin_t^4*cos(phi)^2*sin(phi)^2))^0.5;
                    kf=sqrt(ef*2*m0*q/abs(A)/(1-g))/h;
                    kxf=kf*sin_t*cos(phi);
                    kyf=kf*sin_t*sin(phi);
                    kzf=kf*cos_t;
                    valleyf=3;
                else
                    valleyf=4;
                end
            end            
        elseif r1 < scatGaAs_hole(9,ie,valley)%non-polar, to l, aorbsion
            r=rand();
            flag=0;
            phi=2*pi*rand();
            cos_t=1-2*rand();
            sin_t=sqrt(1-cos_t*cos_t);
            g=((B/A)^2+(C/A)^2*(sin_t*sin_t*cos_t*cos_t+sin_t^4*cos(phi)^2*sin(phi)^2))^0.5;
            fmax=1/(a*(1+B/A))^1.5;
            f=1/(a*(1+g))^1.5;
            i=1;
            while r>f/fmax
                phi=2*pi*rand();
                cos_t=1-2*rand();
                sin_t=sqrt(1-cos_t*cos_t);
                g=((B/A)^2+(C/A)^2*(sin_t*sin_t*cos_t*cos_t+sin_t^4*cos(phi)^2*sin(phi)^2))^0.5;
                f=1/(a*(1+g))^1.5;
                i=i+1;
                if i>2
                    flag=1;
                    break;
                end
            end
            if flag==0
                ef=ei+hw0;
%                 g=((B/A)^2+(C/A)^2*(sin_t^2*cos_t^2+sin_t^4*cos(phi)^2*sin(phi)^2))^0.5;
                kf=sqrt(ef*2*m0*q/abs(A)/(1+g))/h;
                kxf=kf*sin_t*cos(phi);
                kyf=kf*sin_t*sin(phi);
                kzf=kf*cos_t;
                valleyf=4;  
            else
                valleyf=4;  
            end
        elseif r1 < scatGaAs_hole(10,ie,valley)%non-polar, to l, emission
            ef=ei-hw0;
            if ef<0
                valleyf=4;
            else
                r=rand();
                flag=0;
                phi=2*pi*rand();
                cos_t=1-2*rand();
                sin_t=sqrt(1-cos_t*cos_t);
                g=((B/A)^2+(C/A)^2*(sin_t*sin_t*cos_t*cos_t+sin_t^4*cos(phi)^2*sin(phi)^2))^0.5;
%                 a=h*h*abs(A)/(2*m0);
                fmax=1/(a*(1+B/A))^1.5;
                f=1/(a*(1+g))^1.5;
                i=1;
                while r>f/fmax
                    phi=2*pi*rand();
                    cos_t=1-2*rand();
                    sin_t=sqrt(1-cos_t*cos_t);
                    g=((B/A)^2+(C/A)^2*(sin_t*sin_t*cos_t*cos_t+sin_t^4*cos(phi)^2*sin(phi)^2))^0.5;
                    f=1/(a*(1+g))^1.5;
                    i=i+1;
                    if i>2
                        flag=1;
                        break;
                    end
                end
                if flag==0
%                     g=((B/A)^2+(C/A)^2*(sin_t^2*cos_t^2+sin_t^4*cos(phi)^2*sin(phi)^2))^0.5;
                    kf=sqrt(ef*2*m0*q/abs(A)/(1+g))/h;
                    kxf=kf*sin_t*cos(phi);
                    kyf=kf*sin_t*sin(phi);
                    kzf=kf*cos_t;
                    valleyf=4;
                else
                    valleyf=4;
                end            
            end
        elseif r1 < scatGaAs_hole(11,ie,valley)
            r=rand();
            flag=0;
            phi=2*pi*rand();
            cos_t=1-2*rand();
            sin_t=sqrt(1-cos_t*cos_t);
            g=((B/A)^2+(C/A)^2*(sin_t*sin_t*cos_t*cos_t+sin_t^4*cos(phi)^2*sin(phi)^2))^0.5;
%             a=h*h*abs(A)/(2*m0);
            fmax=1/(a*(1+B/A))^1.5;
            f=1/(a*(1+g))^1.5;
            i=1;
            while r>f/fmax
                phi=2*pi*rand();
                cos_t=1-2*rand();
                sin_t=sqrt(1-cos_t*cos_t);
                g=((B/A)^2+(C/A)^2*(sin_t*sin_t*cos_t*cos_t+sin_t^4*cos(phi)^2*sin(phi)^2))^0.5;
                f=1/(a*(1+g))^1.5;
                i=i+1;
                if i>2
                    flag=1;
                    break;
                end
            end
            if flag==0
                ef=ei;
%                 g=((B/A)^2+(C/A)^2*(sin_t^2*cos_t^2+sin_t^4*cos(phi)^2*sin(phi)^2))^0.5;
                kf=sqrt(ef*2*m0*q/abs(A)/(1+g))/h;
                kxf=kf*sin_t*cos(phi);
                kyf=kf*sin_t*sin(phi);
                kzf=kf*cos_t;
                valleyf=4;   
            else
                valleyf=4;   
            end
        elseif r1 < scatGaAs_hole(12,ie,valley)
            r=rand();
            flag=0;
            phi=2*pi*rand();
            cos_t=1-2*rand();
            sin_t=sqrt(1-cos_t*cos_t);
            g=((B/A)^2+(C/A)^2*(sin_t*sin_t*cos_t*cos_t+sin_t^4*cos(phi)^2*sin(phi)^2))^0.5;
%             a=h*h*abs(A)/(2*m0);
            fmax=1/(a*(1+B/A))^1.5;
            f=1/(a*(1+g))^1.5;
            i=1;
            while r>f/fmax
                phi=2*pi*rand();
                cos_t=1-2*rand();
                sin_t=sqrt(1-cos_t*cos_t);
                g=((B/A)^2+(C/A)^2*(sin_t*sin_t*cos_t*cos_t+sin_t^4*cos(phi)^2*sin(phi)^2))^0.5;
                f=1/(a*(1+g))^1.5;
                i=i+1;
                if i>2
                    flag=1;
                    break;
                end
            end
            if flag==0
                ef=ei;
%                 g=((B/A)^2+(C/A)^2*(sin_t^2*cos_t^2+sin_t^4*cos(phi)^2*sin(phi)^2))^0.5;
                kf=sqrt(ef*2*m0*q/abs(A)/(1-g))/h;
                kxf=kf*sin_t*cos(phi);
                kyf=kf*sin_t*sin(phi);
                kzf=kf*cos_t;
                valleyf=3;   
            else 
                valleyf=4;
            end
        end
    end
end            

particle(1,1)=kxf;
particle(1,2)=kyf;
particle(1,3)=kzf;