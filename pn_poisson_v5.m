function [fx,fy,phi,k]=pn_poisson_v5(dx,dy,nx1,ny1,eps_stat,q,p_charge,n_charge,bg_charge,phi,Vp,Vn)

k=0;
delta_phi=0.01;
delta_phi_max=3e-6;
net_charge=bg_charge-n_charge+p_charge;

%for k=0:POISSON_ITER_MAX
while abs(delta_phi)>3e-6
    for j=1:ny1
        phi(j,1)=Vp;
        phi(j,nx1)=Vn;
    end
    %for top and bottom
    for i=2:nx1-1
        phi(1,i)=phi(2,i);
        phi(ny1,i)=phi(ny1-1,i);
    end
    
    for j=1:ny1
        for i=1:nx1
            phi_temp(j,i)=phi(j,i);
        end
    end   
    
    for j=2:ny1-1
        for i=2:nx1-1            
            phi(j,i)=0.25*(phi_temp(j,i+1)+phi(j,i-1)+phi(j-1,i)+phi_temp(j+1,i)+dx*dx*(net_charge(j,i))*q/eps_stat);
            delta_phi_temp=phi(j,i)-phi_temp(j,i);
            if abs(delta_phi_temp)>abs(delta_phi_max)
                delta_phi_max=delta_phi_temp;
            end
        end
    end
    delta_phi=delta_phi_max;
    delta_phi_max=3e-6;
    k=k+1;
end

%after getting the final phi, check the boundary again
%for left and right
for j=1:ny1
    phi(j,1)=Vp;
    phi(j,nx1)=Vn;
end
%for top and bottom
for i=2:nx1-1
    phi(1,i)=phi(2,i);
    phi(ny1,i)=phi(ny1-1,i);
end

fx(ny1,nx1)=0;
for i=1:ny1
    for j=2:(nx1-1)
        fx(i,j)=((phi(i,j-1)-phi(i,j+1))/dx)/2;%924
        %because E=-d(phi)/dx!!
    end
end
%for left and right boundary
fx(:,1)=fx(:,2);
fx(:,end)=fx(:,end-1);

     
fy(ny1,nx1)=0;       
for j=1:nx1
    for i=2:(ny1-1)
        fy(i,j)=(phi(i-1,j)-phi(i+1,j))/(2*dy);        
    end
end
%for top and bottom boundary
fy(1,:)=fy(2,:);
fy(end,:)=fy(end-1,:);


    
            