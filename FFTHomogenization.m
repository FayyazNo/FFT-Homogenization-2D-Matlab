clc
clear all


%Elastic Properties of the Inclusion
Em=300.9e9;
vm=.35;

%Elastic Properties of the Inclusion
Ei=60e7;
vi=0.23;

mum=Em/(1+vm)/2;
lam=vm*Em/(1+vm)/(1-2*vm);

mui=Ei/(1+vi)/2;
lai=vi*Ei/(1+vi)/(1-2*vi);



RVE=~(im2bw((imread('1.png'))));




x=size(RVE,1)
y=size(RVE,2)

% Stiffness Tensors of the Inclusion and Matrix
Ci=[lai+2*mui, lai,       0;
    lai,       lai+2*mui, 0;
    0,         0,         mui];

Cm=[lam+2*mum, lam,       0;
    lam,       lam+2*mum, 0;
    0,         0,         mum];

%Refrence Medium Lamme Coanstants
mu0=10*(mum+mui)/2;
la0=10*(lam+lai)/2;

%--------------------------------------------------
ep   =zeros( x, y, 3);
sigh =zeros( x, y, 3);
eph_o=zeros( x, y, 3);
eph_n=zeros( x, y, 3);
tau  =zeros( x, y, 3);
tauH =zeros( x, y, 3);


E=[0.0005, 0.0005, 0];

ep(:,:,1)=E(1);
ep(:,:,2)=E(2);
ep(:,:,3)=E(3);


 

sig=multpCx(Ci, Cm, RVE, ep);

gam=zeros(3,3,x,y);

% Green Operator stores in (3,3,Nx,Nx) 4-D array
%Gama==========================================================

s1=0;
s2=0;
s3=0;
s4=0;
s5=0;
s6=0;



c1=1/4/mu0;
c2=(la0+mu0)/mu0/(mu0+2*la0);
for i=1:x
     for j=1:y
     if(i<=(x-1)/2+1)  
         i1=i-1;
     else
         i1=-(x-i)-1;
     end
     
     if(j<=(y-1)/2+1)  
         j1=j-1;
     else
         j1=-(y-j)-1;
     end
        
zet=[i1, j1];

gam(1,1,i,j)=(c1/norm(zet)^2)*(4*zet(1)^2)                     -(c2/norm(zet)^4)*(zet(1)^4);             %gam1111
gam(1,2,i,j)=0                                                 -(c2/norm(zet)^4)*(zet(1)^2*zet(2)^2);    %gam1122
gam(1,3,i,j)=((c1/norm(zet)^2)*(2*zet(1)*zet(2))               -(c2/norm(zet)^4)*(zet(1)^3*zet(2)^1)); %gam1112
gam(2,2,i,j)=((c1/norm(zet)^2)*(4*zet(2)^2)                    -(c2/norm(zet)^4)*(zet(2)^4));            %gam2222
gam(2,3,i,j)=((c1/norm(zet)^2)*(zet(1)*zet(2)*2)               -(c2/norm(zet)^4)*(zet(2)^3*zet(1)^1)); %gam2212
gam(3,3,i,j)=((c1/norm(zet)^2)*(zet(1)^2+zet(2)^2)             -(c2/norm(zet)^4)*(zet(1)^2*zet(2)^2)); %gam1212 
gam(2,1,i,j)= gam(1,2,i,j); 
gam(3,1,i,j)= gam(1,3,i,j);
gam(3,2,i,j)= gam(2,3,i,j);
 
if(i==1 && j==1)
gam(1,1,i,j)=0;
gam(1,2,i,j)=0;
gam(1,3,i,j)=0;
gam(2,2,i,j)=0;
gam(2,3,i,j)=0;
gam(3,3,i,j)=0;
gam(2,1,i,j)= gam(1,2,i,j); 
gam(3,1,i,j)= gam(1,3,i,j);
gam(3,2,i,j)= gam(2,3,i,j);
end   

s1=s1+gam(1,1,i,j);
s2=s2+gam(1,2,i,j);
s3=s3+gam(1,3,i,j);
s4=s4+gam(2,2,i,j);
s5=s5+gam(2,3,i,j);
s6=s6+gam(3,3,i,j);

    end
end



% Iteratve Solution. Here Assumed to be 300 Itteration

for t=1:30

sigh(:,:,1)=fft2(sig(:,:,1));
sigh(:,:,2)=fft2(sig(:,:,2));
sigh(:,:,3)=fft2(sig(:,:,3));


for i=1:x
    for j=1:y
        
     eph_n(i,j,:)= reshape(eph_o(i,j,:),3,1) -gam(:,:,i,j)*reshape(sigh(i,j,:),3,1);   

    end
end


eph_n(1,1,:)=[E(1), E(2), E(3)]*x*y;

ep(:,:,1)=ifft2(eph_n(:,:,1));
ep(:,:,2)=ifft2(eph_n(:,:,2));
ep(:,:,3)=ifft2(eph_n(:,:,3));

eph_o=eph_n;

sig=multpCx(Ci, Cm, RVE, ep);

sum1=0;


% Calculation of the Error
for i=1:x
    for j=1:y
     if(i<=(x-1)/2+1)  
         i1=i-1;
     else
         i1=-(x-i)-1;
     end
     
     if(j<=(y-1)/2+1)  
         j1=j-1;
     else
         j1=-(y-j)-1;
     end
        
     zet=[i1, j1];

       sum1=sum1+norm([sigh(i,j,1)*zet(1), sigh(i,j,2)*zet(2), sigh(i,j,3)*(zet(1)+zet(2))])^2;


    end
end

err(t)=norm(sum1)^0.5/norm([sigh(1,1,1), sigh(1,1,2), sigh(1,1,3) ])/x/y;





'ep---------------------'
reshape(ep(1,1,:),3,1)

reshape(sig(1,1,:),3,1)
'----------------------'
end



%Plots===============================================================
a1=real(ep(:,:,1));
a2=real(ep(:,:,2));
a3=real(ep(:,:,3));

figure(10)
ss1=pcolor(a1);
set(ss1, 'EdgeColor', 'none');
title('\epsilon_1_1','fontSize',16)
colorbar();

figure(20)
ss2=pcolor(a2);
set(ss2, 'EdgeColor', 'none');
title('\epsilon_2_2','fontSize',16)
colorbar();

figure(30)
ss3=pcolor(a3);
set(ss3, 'EdgeColor', 'none');
title('\epsilon_1_2','fontSize',16)
colorbar();


b1=real(sig(:,:,1));
b2=real(sig(:,:,2));
b3=real(sig(:,:,3));

figure(40);
ss4=pcolor(b1);
set(ss4, 'EdgeColor', 'none');
title('\sigma_1_1','fontSize',16)
colorbar();

figure(50);
ss5=pcolor(b2);
set(ss5, 'EdgeColor', 'none');
title('\sigma_2_2','fontSize',16)
colorbar();

figure(60)
ss6=pcolor(b3);
set(ss6, 'EdgeColor', 'none');
title('\sigma_1_2','fontSize',16)
colorbar();


figure(70)
plot(err)
xlabel('Itteration') % x-axis label
ylabel('Error') % y-axis label


