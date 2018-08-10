% function [particles,valley,s_enter,s_exit,s_real,d_enter,d_exit,d_real]=gaas_mesfet_renew_v6(particles,valley,Ls,Ld,Ltot,dx,dy,nx1,ny1,max_particles,bk,T,q,h,alpha,eM,GmG,tdt,s_pts,d_pts,bgc_tr,ppsp,s_enter,s_exit,s_real,d_enter,d_exit,d_real,ti,source_temp,drain_temp,s_icpg)
function [particles,valley,p_added,n_added,number]=pn_renew_v6(particles,valley,Ttot,dx,dy,nx1,ny1,max_particles,p_icpg,n_icpg,bk,T,q,h,alpha,eM,emR,Gm,tdt,left_pts,right_pts,Ltot,A,B,C,ti,number,hd)

xmax=(nx1-1)*dx;
ymax=(ny1-1)*dy;

%--------------how many particles are under simulation--------------------
temp1=valley(:,1)==1;
temp2=valley(:,1)==2;
temp3=valley(:,1)==3;
temp4=valley(:,1)==4;
temp1=temp1.*(1:max_particles).';
temp2=temp2.*(1:max_particles).';
temp3=temp3.*(1:max_particles).';
temp4=temp4.*(1:max_particles).';
temp1=temp1(temp1~=0);
temp2=temp2(temp2~=0);
temp3=temp3(temp3~=0);
temp4=temp4(temp4~=0);
number(ti,1)=length(temp1);
number(ti,2)=length(temp2);
number(ti,3)=length(temp3);
number(ti,4)=length(temp4);

%--------------positive contact Calculations-----------------------------------------
%Find Particles In Half Cells at positive contact
temp1=particles(:,5)<0.5*dx;
temp2=valley(:,1)==3;
temp3=valley(:,1)==4; %Only positive particles

temp4=temp1&temp2;
temp5=temp1&temp3;
temp6=valley(:,1)~=9;

ph=temp4|temp5;%only add and delte holes or both holes and electron?
ph=ph.*(1:max_particles).';
ph=ph(ph~=0);

%temp1=particles(:,5)<0.5*dx;
temp2=valley(:,1)==1;
temp3=valley(:,1)==2; %Only positive particles

temp4=temp1&temp2;
temp5=temp1&temp3;
%temp6=valley(:,1)~=9;

pei=temp4|temp5;%only add and delte holes or both holes and electron?
pei=pei.*(1:max_particles).';
pei=pei(pei~=0);

%Count particles belonging to each grid point...0.5*dx to left or right of
%each grid location (except first and last grid pts)
p_cpg(length(left_pts),1)=0;
for i=1:length(ph)
    index=ph(i);
    y=particles(index,6);
    
    %If particle in first half cell...belongs to grid pt 1
    if y < 0.5*dy
        p_cpg(1,1)=p_cpg(1,1)+1;
    %If particle in last half cell...belongs to last grid pt
    elseif y >= Ttot-0.5*dy
        p_cpg(end,1)=p_cpg(end,1)+1;
    else
        yi=round(particles(index,6)/dy)+1;
        p_cpg(yi,1)=p_cpg(yi,1)+1;
    end
end
for i=1:length(pei)
    index=pei(i);
    y=particles(index,6);
    
    %If particle in first half cell...belongs to grid pt 1
    if y < 0.5*dy
        p_cpg(1,1)=p_cpg(1,1)-1;
    %If particle in last half cell...belongs to last grid pt
    elseif y >= Ttot-0.5*dy
        p_cpg(end,1)=p_cpg(end,1)-1;
    else
        yi=round(particles(index,6)/dy)+1;
        p_cpg(yi,1)=p_cpg(yi,1)-1;
    end
end
p_delta=p_cpg-p_icpg;%hd means highly doped factor
p_added=-1*sum(p_delta);%for p area, it will add positive particles, how many holes are added

%Indices of Particles to Be Deleted
v9i=~temp6;
v9i=v9i.*(1:max_particles).';
v9i=v9i(v9i~=0);

v9_count=1;
for i=1:length(p_delta)
    p_dif=abs(p_delta(i));
    %Need to Add Particles 
    if p_delta(i) < 0
        for j=1:p_dif
            index=v9i(v9_count);
            
            iv=4;%only generate light holes
            ei=-(bk*T/q)*log(rand())*1.5;             
            cos_t=1-2*rand();
            sin_t=sqrt(1-cos_t*cos_t);
            phi=2*pi*rand;
            g=((B/A)^2+(C/A)^2*(sin_t^2*cos_t^2+sin_t^4*cos(phi)^2*sin(phi)^2))^0.5;
            ki=sqrt(2*emR*ei*q/(abs(A)*h*h*(1+g)));%start from light
            kx=ki*sin_t*cos(phi);
            ky=ki*sin_t*sin(phi);
            kz=ki*cos_t;
            
            %Want to inject particles so that they are directed into device
            if kx < 0
                kx = -kx;
            end
            
            particles(index,1)=kx;
            particles(index,2)=ky;
            particles(index,3)=kz;
            particles(index,4)=tdt-log(rand())/Gm(iv);
            
            %Only distribute in first half cell for first grid pt
            if i==1
                particles(index,6)=0.5*dy*rand();
            elseif i==length(p_delta)
                particles(index,6)=Ttot-0.5*dy*rand();
            else
                particles(index,6)=(i-1.5)*dy+dy*rand();
            end
            particles(index,5)=0.5*dx*rand();
            valley(index,1)=iv;
            
            v9_count=v9_count+1;
        end
        
    %Need to Delete Particles
    elseif p_delta(i) > 0
        %Get indices of particles belonging to current grid pt
        %If first grid pt, just get particles within first half cell
        if i==1
            psi_temp=particles(ph,6)<0.5*dy;
            psi_temp=psi_temp.*(1:length(ph)).';
            psi_temp=psi_temp(psi_temp~=0);
        %If last grid pt, just get particles within last half cell
        elseif i==length(p_delta)
            psi_temp=particles(ph,6)>(Ttot-0.5*dy);
            psi_temp=psi_temp.*(1:length(ph)).';
            psi_temp=psi_temp(psi_temp~=0);
        else
            psi_temp1=particles(ph,6)>((i-1.5)*dx);
            psi_temp2=particles(ph,6)<((i-0.5)*dx);
            psi_temp=psi_temp1&psi_temp2;
            psi_temp=psi_temp.*(1:length(ph)).';
            psi_temp=psi_temp(psi_temp~=0); 
        end 
        
        %Go through and randomly delete particles
        for j=1:p_dif
            temp_index=1+floor(length(psi_temp)*rand());
            index=psi_temp(temp_index);
            index2=ph(index);
            valley(index2)=9;
            psi_temp(temp_index)=[];
        end
    end
end


%------------negative  Calculations--------------------------------------------
%Find Particles In Half Cells Under negative Grid Pts
temp1=particles(:,5)>Ltot-0.5*dx;
temp2=valley(:,1)==1;
temp3=valley(:,1)==2;
temp4=temp1&temp2;
temp5=temp1&temp3;
temp6=valley(:,1)~=9;

nei=temp4|temp5;
nei=nei.*(1:max_particles).';
nei=nei(nei~=0);

%temp1=particles(:,5)>Ltot-0.5*dx;
temp2=valley(:,1)==3;
temp3=valley(:,1)==4;
temp4=temp1&temp2;
temp5=temp1&temp3;
temp6=valley(:,1)~=9;

npi=temp4|temp5;
npi=npi.*(1:max_particles).';
npi=npi(npi~=0);

%d_cpg: drain count per grid
n_cpg(length(right_pts),1)=0;
for i=1:length(nei)
    index=nei(i);
    y=particles(index,6);
    
    %If particle in first half cell...belongs to grid pt 1
    if y < 0.5*dy
        n_cpg(1,1)=n_cpg(1,1)+1;
    %If particle in last half cell...belongs to last grid pt
    elseif y >= Ttot-0.5*dx
        n_cpg(end,1)=n_cpg(end,1)+1;
    else
        yi=round(particles(index,6)/dy)+1;
        n_cpg(yi,1)=n_cpg(yi,1)+1;
    end
end
for i=1:length(npi)
    index=npi(i);
    y=particles(index,6);
    
    %If particle in first half cell...belongs to grid pt 1
    if y < 0.5*dy
        n_cpg(1,1)=n_cpg(1,1)-1;
    %If particle in last half cell...belongs to last grid pt
    elseif y >= Ttot-0.5*dx
        n_cpg(end,1)=n_cpg(end,1)-1;
    else
        yi=round(particles(index,6)/dy)+1;
        n_cpg(yi,1)=n_cpg(yi,1)-1;
    end
end
n_delta=n_cpg-n_icpg;
n_added=-1*sum(n_delta);%all electrons

%Indices of Particles to Be Deleted
v9i=~temp6;
v9i=v9i.*(1:max_particles).';
v9i=v9i(v9i~=0);

v9_count=1;
for i=1:length(n_delta)
    n_dif=abs(n_delta(i));
    %Need to Add Particles 
    if n_delta(i) < 0
        for j=1:n_dif
            index=v9i(v9_count);
            
            iv=1;%only inject gamma electron
            ei=-(bk*T/q)*log(rand())*1.5;
            ki=(sqrt(2*eM(iv))/h)*sqrt(ei*(1+alpha(iv)*ei))*sqrt(q);
            
            cos_t=1-2*rand();
            sin_t=sqrt(1-cos_t*cos_t);
            phi=2*pi*rand;
            kx=ki*sin_t*cos(phi);
            ky=ki*sin_t*sin(phi);
            kz=ki*cos_t;
            
            %Want to inject particles so that they are directed into device
            if kx > 0
                kx = -kx;
            end
            
            particles(index,1)=kx;
            particles(index,2)=ky;
            particles(index,3)=kz;
            particles(index,4)=tdt-log(rand())/Gm(iv);
            
            %Only distribute in first half cell for first grid pt
            %-----------------------
            if i==1
                particles(index,6)=0.5*dy*rand();
            elseif i==length(n_delta)
                particles(index,6)=Ttot-0.5*dy*rand();
            else
                particles(index,6)=(i-1.5)*dy+dy*rand();
            end
            particles(index,5)=Ltot-0.5*dx*rand();
            valley(index,1)=iv;
            
            v9_count=v9_count+1;
        end
        
    %Need to Delete Particles
    elseif n_delta(i) > 0
        %Get indices of particles belonging to current grid pt
        %If first grid pt, just get particles within first half cell
        if i==1
            pdi_temp=particles(nei,6)<0.5*dy;
            pdi_temp=pdi_temp.*(1:length(nei)).';
            pdi_temp=pdi_temp(pdi_temp~=0);
        %If last grid pt, just get particles within last half cell
        elseif i==length(n_delta)
            pdi_temp=particles(nei,5)>(Ttot-0.5*dy);
            pdi_temp=pdi_temp.*(1:length(nei)).';
            pdi_temp=pdi_temp(pdi_temp~=0);
        else
            pdi_temp1=particles(nei,6)>((i-1.5)*dx);
            pdi_temp2=particles(nei,6)<((i-0.5)*dx);
            pdi_temp=pdi_temp1&pdi_temp2;
            pdi_temp=pdi_temp.*(1:length(nei)).';
            pdi_temp=pdi_temp(pdi_temp~=0); 
        end 
        
        %Go through and randomly delete particles
        for j=1:n_dif
            temp_index=ceil(length(pdi_temp)*rand());
            index=pdi_temp(temp_index);
            index2=nei(index);
            valley(index2)=9;
            pdi_temp(temp_index)=[];
        end
    end
end





    
