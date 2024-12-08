global gs5 wt5 gsx gsy wt  mu lambda tf theta  a b ko  gps1 gps2  Wr W Wrr xRight
global gpL gpR gpB gpT  gcL1 gcR1 gcB1 gcT1  gcL2 gcR2 gcB2 gcT2 gpM tfx tfy aai  myeps
% mu=1/512; lambda=1/65500; 
mu=1/400;  lambda=1/40000;   
xRight=1;   myeps=1.0e-5;   ko=3;
gk=5;  theta=0.85; 
gs5=[-0.9061798459, 0.9061798459, -0.5384693101, 0.5384693101, 0.0];
wt5=[0.2369268851, 0.2369268851, 0.4786286705,  0.4786286705, 0.5688888889]';
% a=8.0*0.85/27.0;   b=1.0;   
% Wrr=@(r) a*b^2./(r.*(b - r).^2)-2;
% Wr=@(r) a*(log(r./(b-r))+b./(b-r))-2*r;
% W=@(r) a*r.*log(r./(b - r)) - r.^2;
a=3.0;   b=3.0;  R=8.0/3.0;   
Wrr=@(r) -(2*a*b^2*r - R*theta*b^2 - 4*a*b*r.^2 + 2*a*r.^3)./(r.*(b - r).^2);
Wr=@(r) R*theta*log(r./(b - r)) - 2*a*r + (R*b*theta)./(b - r);
W=@(r) R*theta*r.*log(r./(b-r))-a*r.^2;

id=0;   wt=zeros(gk*gk,1); gsx=zeros(1,gk*gk); gsy=zeros(1,gk*gk);
for i=1:gk
    for j=1:gk
        id=id+1;
        wt(id)=wt5(i)*wt5(j);   gsx(id)=gs5(i);    gsy(id)=gs5(j);
    end
end  
gpR=ones(ko,gk);    gpR(2,:)=-1.0;     gpR(3,:)=gs5;
gpL=ones(ko,gk);    gpL(3,:)=gs5;
gpT=ones(ko,gk);    gpT(2,:)=gs5;      gpT(3,:)=-1.0; 
gpB=ones(ko,gk);    gpB(2,:)=gs5;
gcR1=gpR;   gcL1=gpL;   gcT1=gpT;   gcB1=gpB;
gcR2=gpR;   gcL2=gpL;   gcT2=gpT;   gcB2=gpB;
gps1=0.5*(gs5-1.0);  gps2=0.5*(gs5+1.0);
gcR1(3,:)=gps1;    gcL1(3,:)=gps1;  gcT1(2,:)=gps1;    gcB1(2,:)=gps1;
gcR2(3,:)=gps2;    gcL2(3,:)=gps2;  gcT2(2,:)=gps2;    gcB2(2,:)=gps2;
oneg=ones(1,gk*gk); zerog=zeros(1,gk*gk); 
tfx=[zerog; oneg; zerog];  tfy=[zerog; zerog; oneg];


ko=3;   gk=5;
gs5=[-0.9061798459, 0.9061798459, -0.5384693101, 0.5384693101, 0.0];
wt5=[0.2369268851, 0.2369268851, 0.4786286705,  0.4786286705, 0.5688888889]';

id=0;   wt=zeros(gk*gk,1); gsx=zeros(1,gk*gk); gsy=zeros(1,gk*gk);
for i=1:gk
    for j=1:gk
        id=id+1;
        wt(id)=wt5(i)*wt5(j);   gsx(id)=gs5(i);    gsy(id)=gs5(j);
    end
end  

gpM=[ones(1,gk*gk); gsx;  gsy];

c2f1=zeros(ko,ko);  c2f2=c2f1;  c2f3=c2f1;  c2f4=c2f1;  tf=gpM;
xi1=(gsx-1)/2;  xi2=(gsx+1)/2;  yi1=(gsy-1)/2;  yi2=(gsy+1)/2;
tf1=[ones(1,gk*gk); xi1;  yi1]; tf2=[ones(1,gk*gk); xi2;  yi1];
tf3=[ones(1,gk*gk); xi2;  yi2]; tf4=[ones(1,gk*gk); xi1;  yi2];
for i=1:ko
    for j=1:ko
        c2f1(i,j)=tf1(j,:).*tf(i,:)*wt; c2f2(i,j)=tf2(j,:).*tf(i,:)*wt;
        c2f3(i,j)=tf3(j,:).*tf(i,:)*wt; c2f4(i,j)=tf4(j,:).*tf(i,:)*wt;
    end
end
aai=diag([1 3 3]);
c2f1=aai*c2f1/4;    c2f2=aai*c2f2/4;    c2f3=aai*c2f3/4;    c2f4=aai*c2f4/4;
