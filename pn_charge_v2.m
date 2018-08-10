function [charge_p,charge_n]=pn_charge_v2(particles,valley,nx1,ny1,dx,dy,max_particles,cpsp)

charge_p(ny1,nx1)=0;
charge_n(ny1,nx1)=0;
temp_n=0;

for n=1:max_particles
    if valley(n)~=9
        x=particles(n,5)/dx;
        y=particles(n,6)/dy;

        i=min(floor(y)+1,ny1-1);
        j=min(floor(x)+1,nx1-1);       
        
        i=max(i,1);
        j=max(j,1);

        yb=i-1;
        xb=j-1;

        y1=1-(y-yb);
        x1=1-(x-xb);
        if valley(n)==1||valley(n)==2
            charge_n(i,j)=charge_n(i,j)+x1*y1;
            charge_n(i+1,j)=charge_n(i+1,j)+x1*(1-y1);%final check, correct 924
            charge_n(i,j+1)=charge_n(i,j+1)+y1*(1-x1);
            charge_n(i+1,j+1)=charge_n(i+1,j+1)+(1-x1)*(1-y1);
            temp_n=temp_n+1;%count how many particles under simulation
        elseif valley(n)==3||valley(n)==4
%             jj='n=';%914
%             fprintf('%s %f\n',jj,n);%914
            charge_p(i,j)=charge_p(i,j)+x1*y1;
            charge_p(i+1,j)=charge_p(i+1,j)+x1*(1-y1);
            charge_p(i,j+1)=charge_p(i,j+1)+y1*(1-x1);
            charge_p(i+1,j+1)=charge_p(i+1,j+1)+(1-x1)*(1-y1);
            temp_n=temp_n+1;
        end
    end
end

for i=1:ny1
    for j=1:nx1
        charge_p(i,j)=charge_p(i,j)*cpsp/dx/dy;
        charge_n(i,j)=charge_n(i,j)*cpsp/dx/dy;
        if i==1 || i==ny1
            charge_p(i,j)=charge_p(i,j)*2;
            charge_n(i,j)=charge_n(i,j)*2;
        end
        
        if j==1 || j==nx1
            charge_p(i,j)=charge_p(i,j)*2;
            charge_n(i,j)=charge_n(i,j)*2;
        end
    end
end

charge_p(ny1,nx1)=charge_p(ny1,nx1-1);%928,from JD,top-right
charge_p(ny1,1)=charge_p(ny1-1,1);%top-left
charge_n(ny1,nx1)=charge_n(ny1,nx1-1);%928,from JD,top-right
charge_n(ny1,1)=charge_n(ny1-1,1);%top-left




        