load('ig=50id=0.mat');

x_charge(max_particles,1)=0;
result(max_particles,1)=0;

for n=1:max_particles
    iv=valley(n,1);
    if iv==1||iv==2
        kx=particles(n,1);
        ky=particles(n,2);
        kz=particles(n,3);

        skx=kx*kx;
        sky=ky*ky;
        skz=kz*kz;
        sk=abs(skx+sky+skz);
        ki=sqrt(sk);

        cos_theta=kz/ki;
        sin_theta=sqrt(1-cos_theta^2);
        cos_phi=kx/(ki*sin_theta);
        sin_phi=ky/(ki*sin_theta);
        g1=((B/A)^2+(C/A)^2*(sin_theta^2*cos_theta^2+sin_theta^4*cos_phi^2*sin_phi^2))^0.5;

        gk=(h*h/(2*eM(iv)))*sk*(1/q);
        if iv==1||iv==2
            ei=(sqrt(1+4*alpha(iv)*gk)-1)/(2*alpha(iv));
            particles(n,7)=ei;
        elseif iv==3
            ei=(h*h*abs(A)/(2*emR))*sk*(1/q)*(1-g1);            
            particles(n,7)=ei;
        elseif iv==4
            ei=(h*h*abs(A)/(2*emR))*sk*(1/q)*(1+g1);
            particles(n,7)=ei;
        end
        %-----connect position and energy together-----
        result(n,1)=particles(n,5);
        result(n,2)=particles(n,7);
    end
end

i=1;
temp1=result(:,1)~=0;
for n=1:max_particles    
    if temp1(n,1)==1
        ar(i,:)=result(n,:);
        i=i+1;
    end
end

figure
scatter(ar(:,1),ar(:,2));



    
    

