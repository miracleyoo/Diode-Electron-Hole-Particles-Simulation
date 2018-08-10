function [particles,valley,bg_charge,cpsp,N,p_icpg,n_icpg,xmax,ymax]=pn_init_v2(max_particles,dope_type,dx,dy,Ltot,nx1,ny1,cdop,ppc,bk,T,q,h,alpha,eM,Gm,Lp,A,B,C,emR,hd)
%dope_type:  1:Uniform 2:Only S/D 3:Uniform+Contacts 4:S/D+Contacts
%ppsp:       particles per super particle

%Create Matrix/Vector for Particles and Valley Tracking 
particles(max_particles,6)=0;
valley(max_particles,1)=0;
valley(:,1)=9;  %Set all valley indices to 9 initially...will be set properly below
N=0;

p_x=Lp/dx+1;

%Create Doping Profile
bg_charge(ny1,nx1)=0;
if dope_type==1
    for i=1:ny1
        for j=1:nx1
            bg_charge(i,j)=cdop;
            if j<p_x
                bg_charge(i,j)=-1*cdop;
            end
        end
    end
end

%Calculate cpsp (charge per super particle) based on requested ppc (particles
%per cell) for the highest doped cell
cpsp=(cdop*dx*dy)/ppc;

xmax=(nx1-1)*dx;
ymax=(ny1-1)*dy;
n=0;
for i=1:ny1
    for j=1:nx1
        npij=floor(abs(bg_charge(i,j))*dx*dy/cpsp+0.5);
        
        %Cells Around the Outside only have half the area, so they only get
        %half the number of particles
        if i==1 || i==ny1
            npij=npij/2;
        end
        if j==1 || j==nx1
            npij=npij/2*hd;
        end
        
        for m=1:npij
            n=n+1;
            if n > max_particles
                disp('Max Number of Particles Exceeded!');
                break;
            end
            
            if j<p_x%positive area
                iv=4;%start from light holes
            else
                iv=1;
            end
            ei=-(bk*T/q)*log(rand())*1.5;
            if iv==4              
                cos_t=1-2*rand();
                sin_t=sqrt(1-cos_t*cos_t);
                phi=2*pi*rand;
                g=((B/A)^2+(C/A)^2*(sin_t^2*cos_t^2+sin_t^4*cos(phi)^2*sin(phi)^2))^0.5;
                ki=sqrt(2*emR*ei*q/(abs(A)*h*h*(1+g)));%start from heavy
                kx=ki*sin_t*cos(phi);
                ky=ki*sin_t*sin(phi);
                kz=ki*cos_t;                
            else
                ki=(sqrt(2*eM(iv))/h)*sqrt(ei*(1+alpha(iv)*ei))*sqrt(q);
                cos_t=1-2*rand();
                sin_t=sqrt(1-cos_t*cos_t);
                phi=2*pi*rand;
                kx=ki*sin_t*cos(phi);
                ky=ki*sin_t*sin(phi);
                kz=ki*cos_t;
            end            
            
            particles(n,1)=kx;                  %kx
            particles(n,2)=ky;                  %ky
            particles(n,3)=kz;                  %kz   
            particles(n,4)=-log(rand())/Gm(iv);    %ts
            particles(n,5)=dx*(rand()+j-1.5);   %x
            particles(n,6)=dy*(rand()+i-1.5);   %z
            particles(n,7)=0;
            
            %------------------------------------------
            if i==1
                particles(n,6)=dy*0.5*rand();
            end
            
            if j==1
                particles(n,5)=dx*0.5*rand();
            end
            
            if i==ny1
                particles(n,6)=ymax-dy*0.5*rand();
            end
            
            if j==nx1
                particles(n,5)=xmax-dx*0.5*rand();
            end
            
            valley(n,1)=iv;
            N=N+1;
        end
    end
end

%Number of Particles under p contact
temp1=particles(:,5)<0.5*dx;
temp3=valley(:,1)~=9;

psi=temp1&temp3;
p_total=sum(psi); 
psi=psi.*(1:max_particles).';
psi=psi(psi~=0);

%--------------------Added 11/1------------------------------------
%Count particles belonging to each grid point...0.5*dx to left or right of
%each grid location (except first and last grid pts)
p_icpg=zeros(ny1,1);
for i=1:length(psi)
    index=psi(i);
    y=particles(index,6);

    %If particle in first half cell...belongs to grid pt 1
    if y < 0.5*dy
        p_icpg(1,1)=p_icpg(1,1)+1;
    %If particle in last half cell...belongs to last grid pt
    elseif y >= Ltot-0.5*dx
        p_icpg(end,1)=p_icpg(end,1)+1;
    else
        yi=round(particles(index,6)/dy)+1;
        p_icpg(yi,1)=p_icpg(yi,1)+1;
    end
end
%------------------------------------------------------------------ 
%Number of Particles near negative contact
temp1=particles(:,5)>(xmax-0.5*dx);
pdi=temp1&temp3;
n_total=sum(pdi);
pdi=pdi.*(1:max_particles).';
pdi=pdi(pdi~=0);

%--------------------Added 11/1------------------------------------ 
n_icpg=zeros(ny1,1);
for i=1:length(pdi)
    index=pdi(i);
    y=particles(index,6);

    %If particle in first half cell...belongs to grid pt 1
    if y < 0.5*dy
        n_icpg(1,1)=n_icpg(1,1)+1;
    %If particle in last half cell...belongs to last grid pt
    elseif y >= Ltot-0.5*dx
        n_icpg(end,1)=n_icpg(end,1)+1;
    else
        xi=round(particles(index,6)/dy)+1;
        n_icpg(xi,1)=n_icpg(xi,1)+1;
    end
end