function Mat=Jacobx(Cellsit, Mesh, Mesh_fine, dt)
global gs5 wt5 gps1 gps2 wt mu lambda  aai Wrr Wr W tfx tfy myeps 
global gpL gpR gpB gpT  gcL1 gcR1 gcB1 gcT1  gcL2 gcR2 gcB2 gcT2 gpM
ko=3;    nk=10*ko;  gk=5;  xxi=gpM(2,:);   yyi=gpM(3,:);  aav=diag(aai);
%%%%%%%%%%%%%%%%%%%%% Mesh_fine u1,u2,p
fCells=Mesh_fine.fCells;    fCellsn=Mesh_fine.fCellsn;
fEy=Mesh_fine.fEy;   fEyn=Mesh_fine.fEyn;   
fEx=Mesh_fine.fEx;   fExn=Mesh_fine.fExn;             
Cells=Mesh.Cells;    Ey=Mesh.Ey;     Ex=Mesh.Ex; 
nks=fCellsn*nk; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% q_1=rho_x, q_2=rho_y %%rho_t+(q_1)_x+(q_2)_y
ko2=ko*ko;  oneb=ones(1,gk);    
kfy=0;    Ify=zeros(1,12*nks);     Jfy=Ify;     sfy=Ify;
kfye=0;  Ifye=zeros(1,12*nks);   Jfye=Ifye;  sfye=Ifye;
for k=1:fEyn 
    id=fEy(k);   EK=Ey(id);   K1=Cells(EK.K1);    K2=Cells(EK.K2);  
    K1it=Cellsit(EK.K1);    K2it=Cellsit(EK.K2);   h1=K1.Wid;  h2=K2.Wid;
    EK1f=K1.Indx_fine;  EK2f=K2.Indx_fine;  h=min(h1, h2); h1e=h/h1/h1; h2e=h/h2/h2;
    if K1.Lev==K2.Lev
         xib1=gs5;   ULn=K1.U*gpL;   UL=K1it.U*gpL; 
         xib2=gs5;   URn=K2.U*gpR;   UR=K2it.U*gpR; 
    elseif K1.Lev==K2.Lev-1
        if K1.Center(2)>K2.Center(2)
            xib1=gps1;  ULn=K1.U*gcL1; UL=K1it.U*gcL1;   
        else
            xib1=gps2;  ULn=K1.U*gcL2; UL=K1it.U*gcL2;   
        end
        xib2=gs5;   URn=K2.U*gpR;   UR=K2it.U*gpR;  
    elseif K1.Lev==K2.Lev+1
        if K2.Center(2)>K1.Center(2)
            URn=K2.U*gcR1;  xib2=gps1;  UR=K2it.U*gcR1;  
        else
            URn=K2.U*gcR2;  xib2=gps2;  UR=K2it.U*gcR2; 
        end
        xib1=gs5;   ULn=K1.U*gpL;    UL=K1it.U*gpL;   
    else
        disp('Init: wrong conforming')
    end  
    rRn=URn(1,:);  uRn=URn(2,:);  rR=UR(1,:);  uR=UR(2,:);  zR=UR(4,:);
    vRn=URn(3,:);  vLn=ULn(3,:);  vR=UR(3,:);  vL=UL(3,:);  zL=UL(4,:); 
    
    nEK1 =(EK1f-1)*nk;   nEK2 =(EK2f-1)*nk;
    Iy1l=(nEK1+1)*ones(1,ko);   Iy1r=(nEK2+1)*ones(1,ko);   
    Jy1l=nEK1+1:nEK1+ko;    Jy1r=nEK2+1:nEK2+ko;  
%     %%%%%%%%%%%%%%%%%%%%%% r_t+(ru)_x
    s2e=(uR+uRn)/4;  
    s2   =[s2e*wt5, -s2e*wt5, s2e.*xib2*wt5]*0.5;
    s12xi=[s2e.*xib1*wt5, -s2e.*xib1*wt5, s2e.*xib2.*xib1*wt5]*0.5;
    s2xi =[s2e.*xib2*wt5, -s2e.*xib2*wt5, s2e.*xib2.*xib2*wt5]*0.5;
    Ify(kfy*ko2+1:kfy*ko2+ko2)=[Iy1l, Iy1l+1, Iy1l+2];
    Jfy(kfy*ko2+1:kfy*ko2+ko2)=[Jy1r, Jy1r, Jy1r];
    sfy(kfy*ko2+1:kfy*ko2+ko2)=dt*[s2, 3*s2, 3*s12xi]*h1e;  kfy=kfy+1;
    Ify(kfy*ko2+1:kfy*ko2+ko2)=[Iy1r, Iy1r+1, Iy1r+2];
    Jfy(kfy*ko2+1:kfy*ko2+ko2)=[Jy1r, Jy1r, Jy1r];
    sfy(kfy*ko2+1:kfy*ko2+ko2)=dt*[-s2, 3*s2, -3*s2xi]*h2e; kfy=kfy+1;
    s2e=(rR+rRn)/4;  
    s2   =[s2e*wt5, -s2e*wt5, s2e.*xib2*wt5]*0.5;
    s12xi=[s2e.*xib1*wt5, -s2e.*xib1*wt5, s2e.*xib2.*xib1*wt5]*0.5;
    s2xi =[s2e.*xib2*wt5, -s2e.*xib2*wt5, s2e.*xib2.*xib2*wt5]*0.5;
    Ify(kfy*ko2+1:kfy*ko2+ko2)=[Iy1l, Iy1l+1, Iy1l+2];
    Jfy(kfy*ko2+1:kfy*ko2+ko2)=[Jy1r, Jy1r, Jy1r]+ko;
    sfy(kfy*ko2+1:kfy*ko2+ko2)=dt*[s2, 3*s2, 3*s12xi]*h1e;  kfy=kfy+1;
    Ify(kfy*ko2+1:kfy*ko2+ko2)=[Iy1r, Iy1r+1, Iy1r+2];
    Jfy(kfy*ko2+1:kfy*ko2+ko2)=[Jy1r, Jy1r, Jy1r]+ko;
    sfy(kfy*ko2+1:kfy*ko2+ko2)=dt*[-s2, 3*s2, -3*s2xi]*h2e; kfy=kfy+1;
             
%     rZv=(rR+rRn)/2.*((zR-vR.^2/4-vRn.^2/4)-(zL-vL.^2/4-vLn.^2/4));
%     F2f(EK2f,:)=F2f(EK2f,:)+[rZv*wt5,-rZv*wt5, rZv.*xib2*wt5]*0.5*h; 
    s2e=1/2.*((zR-vR.^2/4-vRn.^2/4)-(zL-vL.^2/4-vLn.^2/4));
    s2   =[s2e*wt5, -s2e*wt5, s2e.*xib2*wt5]*0.5;
    s2xi =[s2e.*xib2*wt5, -s2e.*xib2*wt5, s2e.*xib2.*xib2*wt5]*0.5;
    Ify(kfy*ko2+1:kfy*ko2+ko2)=[Iy1r, Iy1r+1, Iy1r+2]+ko;
    Jfy(kfy*ko2+1:kfy*ko2+ko2)=[Jy1r, Jy1r, Jy1r];
    sfy(kfy*ko2+1:kfy*ko2+ko2)=dt*[s2, -s2, s2xi]*h2e; kfy=kfy+1; 
    s2e=(rR+rRn)/2;
    s2   =[s2e*wt5, -s2e*wt5, s2e.*xib2*wt5]*0.5;
    s2xi =[s2e.*xib2*wt5, -s2e.*xib2*wt5, s2e.*xib2.*xib2*wt5]*0.5;
    Ify(kfy*ko2+1:kfy*ko2+ko2)=[Iy1r, Iy1r+1, Iy1r+2]+ko;
    Jfy(kfy*ko2+1:kfy*ko2+ko2)=[Jy1r, Jy1r, Jy1r]+3*ko;
    sfy(kfy*ko2+1:kfy*ko2+ko2)=dt*[s2, -s2, s2xi]*h2e; kfy=kfy+1; 
    s1e=-(rR+rRn)/2;  
    s1  =[s1e*wt5, s1e*wt5, s1e.*xib1*wt5]*0.5;  
    s21xi=[s1e.*xib2*wt5, s1e.*xib2*wt5, s1e.*xib1.*xib2*wt5]*0.5;
    Ify(kfy*ko2+1:kfy*ko2+ko2)=[Iy1r, Iy1r+1, Iy1r+2]+ko;
    Jfy(kfy*ko2+1:kfy*ko2+ko2)=[Jy1l, Jy1l, Jy1l]+3*ko;
    sfy(kfy*ko2+1:kfy*ko2+ko2)=dt*[s1, -s1, s21xi]*h2e; kfy=kfy+1;
    s2e=(rR+rRn)/2.*(-vR/2);
    s2   =[s2e*wt5, -s2e*wt5, s2e.*xib2*wt5]*0.5;
    s2xi =[s2e.*xib2*wt5, -s2e.*xib2*wt5, s2e.*xib2.*xib2*wt5]*0.5;
    Ify(kfy*ko2+1:kfy*ko2+ko2)=[Iy1r, Iy1r+1, Iy1r+2]+ko;
    Jfy(kfy*ko2+1:kfy*ko2+ko2)=[Jy1r, Jy1r, Jy1r]+2*ko;
    sfy(kfy*ko2+1:kfy*ko2+ko2)=dt*[s2, -s2, s2xi]*h2e; kfy=kfy+1;
    s1e=(rR+rRn)/2.*(vL/2);  
    s1  =[s1e*wt5, s1e*wt5, s1e.*xib1*wt5]*0.5;  
    s21xi=[s1e.*xib2*wt5, s1e.*xib2*wt5, s1e.*xib1.*xib2*wt5]*0.5;
    Ify(kfy*ko2+1:kfy*ko2+ko2)=[Iy1r, Iy1r+1, Iy1r+2]+ko;
    Jfy(kfy*ko2+1:kfy*ko2+ko2)=[Jy1l, Jy1l, Jy1l]+2*ko;
    sfy(kfy*ko2+1:kfy*ko2+ko2)=dt*[s1, -s1, s21xi]*h2e; kfy=kfy+1;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     rZu=(rR+rRn)/8.*(uR+uRn).*(vR+vRn-vL-vLn);
%     F3f(EK2f,:)=F3f(EK2f,:)+[rZu*wt5,-rZu*wt5, rZu.*xib2*wt5]*0.5*h; 
    s2e=1/8.*(uR+uRn).*(vR+vRn-vL-vLn);
    s2   =[s2e*wt5, -s2e*wt5, s2e.*xib2*wt5]*0.5;
    s2xi =[s2e.*xib2*wt5, -s2e.*xib2*wt5, s2e.*xib2.*xib2*wt5]*0.5;
    Ify(kfy*ko2+1:kfy*ko2+ko2)=[Iy1r, Iy1r+1, Iy1r+2]+2*ko;
    Jfy(kfy*ko2+1:kfy*ko2+ko2)=[Jy1r, Jy1r, Jy1r];
    sfy(kfy*ko2+1:kfy*ko2+ko2)=dt*[s2, -s2, s2xi]*h2e; kfy=kfy+1;
    s2e=1/8.*(rR+rRn).*(vR+vRn-vL-vLn);
    s2   =[s2e*wt5, -s2e*wt5, s2e.*xib2*wt5]*0.5;
    s2xi =[s2e.*xib2*wt5, -s2e.*xib2*wt5, s2e.*xib2.*xib2*wt5]*0.5;
    Ify(kfy*ko2+1:kfy*ko2+ko2)=[Iy1r, Iy1r+1, Iy1r+2]+2*ko;
    Jfy(kfy*ko2+1:kfy*ko2+ko2)=[Jy1r, Jy1r, Jy1r]+ko;
    sfy(kfy*ko2+1:kfy*ko2+ko2)=dt*[s2, -s2, s2xi]*h2e; kfy=kfy+1;
    s2e=1/8.*(rR+rRn).*(uR+uRn);
    s2   =[s2e*wt5, -s2e*wt5, s2e.*xib2*wt5]*0.5;
    s2xi =[s2e.*xib2*wt5, -s2e.*xib2*wt5, s2e.*xib2.*xib2*wt5]*0.5;
    Ify(kfy*ko2+1:kfy*ko2+ko2)=[Iy1r, Iy1r+1, Iy1r+2]+2*ko;
    Jfy(kfy*ko2+1:kfy*ko2+ko2)=[Jy1r, Jy1r, Jy1r]+2*ko;
    sfy(kfy*ko2+1:kfy*ko2+ko2)=dt*[s2, -s2, s2xi]*h2e; kfy=kfy+1;
    s1e=-s2e;
    s1  =[s1e*wt5, s1e*wt5, s1e.*xib1*wt5]*0.5;  
    s21xi=[s1e.*xib2*wt5, s1e.*xib2*wt5, s1e.*xib1.*xib2*wt5]*0.5;
    Ify(kfy*ko2+1:kfy*ko2+ko2)=[Iy1r, Iy1r+1, Iy1r+2]+2*ko;
    Jfy(kfy*ko2+1:kfy*ko2+ko2)=[Jy1l, Jy1l, Jy1l]+2*ko;
    sfy(kfy*ko2+1:kfy*ko2+ko2)=dt*[s1, -s1, s21xi]*h2e; kfy=kfy+1;
    %%%%%%%%%%%%%%%%%%%%%%%% Z+lambda*q1_x
    s1e=oneb;   
    s1  =[s1e*wt5, s1e*wt5, s1e.*xib1*wt5]*0.5;  
    s1xi=[s1e.*xib1*wt5, s1e.*xib1*wt5, s1e.*xib1.*xib1*wt5]*0.5;
    s21xi=[s1e.*xib2*wt5, s1e.*xib2*wt5, s1e.*xib1.*xib2*wt5]*0.5;
% %     f2x=-mu*(4/3*Ukxi(7,:)-2/3*Ukxi(10,:)); 
    Ify(kfy*ko2+1:kfy*ko2+ko2)=[Iy1l, Iy1l+1, Iy1l+2]+ko;
    Jfy(kfy*ko2+1:kfy*ko2+ko2)=[Jy1l, Jy1l, Jy1l]+6*ko;
    sfy(kfy*ko2+1:kfy*ko2+ko2)=-4/3*mu*dt*[s1, s1, s1xi]*h1e;   kfy=kfy+1;
    Ify(kfy*ko2+1:kfy*ko2+ko2)=[Iy1r, Iy1r+1, Iy1r+2]+ko;
    Jfy(kfy*ko2+1:kfy*ko2+ko2)=[Jy1l, Jy1l, Jy1l]+6*ko;
    sfy(kfy*ko2+1:kfy*ko2+ko2)=-4/3*mu*dt*[-s1, s1,-s21xi]*h2e; kfy=kfy+1;
    
    Ify(kfy*ko2+1:kfy*ko2+ko2)=[Iy1l, Iy1l+1, Iy1l+2]+ko;
    Jfy(kfy*ko2+1:kfy*ko2+ko2)=[Jy1l, Jy1l, Jy1l]+9*ko;
    sfy(kfy*ko2+1:kfy*ko2+ko2)=2/3*mu*dt*[s1, s1, s1xi]*h1e;   kfy=kfy+1;
    Ify(kfy*ko2+1:kfy*ko2+ko2)=[Iy1r, Iy1r+1, Iy1r+2]+ko;
    Jfy(kfy*ko2+1:kfy*ko2+ko2)=[Jy1l, Jy1l, Jy1l]+9*ko;
    sfy(kfy*ko2+1:kfy*ko2+ko2)=2/3*mu*dt*[-s1, s1,-s21xi]*h2e; kfy=kfy+1;
%     f3x=-mu*(Ukxi(8,:)+Ukxi(9,:));
    Ify(kfy*ko2+1:kfy*ko2+ko2)=[Iy1l, Iy1l+1, Iy1l+2]+2*ko;
    Jfy(kfy*ko2+1:kfy*ko2+ko2)=[Jy1l, Jy1l, Jy1l]+7*ko;
    sfy(kfy*ko2+1:kfy*ko2+ko2)=-mu*dt*[s1, s1, s1xi]*h1e;   kfy=kfy+1;
    Ify(kfy*ko2+1:kfy*ko2+ko2)=[Iy1r, Iy1r+1, Iy1r+2]+2*ko;
    Jfy(kfy*ko2+1:kfy*ko2+ko2)=[Jy1l, Jy1l, Jy1l]+7*ko;
    sfy(kfy*ko2+1:kfy*ko2+ko2)=-mu*dt*[-s1, s1,-s21xi]*h2e; kfy=kfy+1;
    
    Ify(kfy*ko2+1:kfy*ko2+ko2)=[Iy1l, Iy1l+1, Iy1l+2]+2*ko;
    Jfy(kfy*ko2+1:kfy*ko2+ko2)=[Jy1l, Jy1l, Jy1l]+8*ko;
    sfy(kfy*ko2+1:kfy*ko2+ko2)=-mu*dt*[s1, s1, s1xi]*h1e;   kfy=kfy+1;
    Ify(kfy*ko2+1:kfy*ko2+ko2)=[Iy1r, Iy1r+1, Iy1r+2]+2*ko;
    Jfy(kfy*ko2+1:kfy*ko2+ko2)=[Jy1l, Jy1l, Jy1l]+8*ko;
    sfy(kfy*ko2+1:kfy*ko2+ko2)=-mu*dt*[-s1, s1,-s21xi]*h2e; kfy=kfy+1;
    
    Ify(kfy*ko2+1:kfy*ko2+ko2)=[Iy1l, Iy1l+1, Iy1l+2]+3*ko;
    Jfy(kfy*ko2+1:kfy*ko2+ko2)=[Jy1l, Jy1l, Jy1l]+4*ko;
    sfy(kfy*ko2+1:kfy*ko2+ko2)=lambda*[s1, 3.0*s1, 3.0*s1xi]*h1e;   kfy=kfy+1;
    Ify(kfy*ko2+1:kfy*ko2+ko2)=[Iy1r, Iy1r+1, Iy1r+2]+3*ko;
    Jfy(kfy*ko2+1:kfy*ko2+ko2)=[Jy1l, Jy1l, Jy1l]+4*ko;
    sfy(kfy*ko2+1:kfy*ko2+ko2)=lambda*[-s1, 3.0*s1,-3.0*s21xi]*h2e; kfy=kfy+1;
%     %%%%%%%%%%%%%%%%%%%% q1-(r^{n+1}+r^n)_x/2
    s2e=oneb;
    s2   =[s2e*wt5, -s2e*wt5, s2e.*xib2*wt5]*0.5;
    s12xi=[s2e.*xib1*wt5, -s2e.*xib1*wt5, s2e.*xib2.*xib1*wt5]*0.5;
    s2xi =[s2e.*xib2*wt5, -s2e.*xib2*wt5, s2e.*xib2.*xib2*wt5]*0.5;
    Ify(kfy*ko2+1:kfy*ko2+ko2)=[Iy1l, Iy1l+1, Iy1l+2]+4*ko;
    Jfy(kfy*ko2+1:kfy*ko2+ko2)=[Jy1r, Jy1r, Jy1r];
    sfy(kfy*ko2+1:kfy*ko2+ko2)=-0.5*[s2, 3*s2, 3*s12xi]*h1e;  kfy=kfy+1;
    Ify(kfy*ko2+1:kfy*ko2+ko2)=[Iy1r, Iy1r+1, Iy1r+2]+4*ko;
    Jfy(kfy*ko2+1:kfy*ko2+ko2)=[Jy1r, Jy1r, Jy1r];
    sfy(kfy*ko2+1:kfy*ko2+ko2)=-0.5*[-s2, 3*s2, -3*s2xi]*h2e; kfy=kfy+1; 
    %%%%%%%%%%%%%%%%%%%% q3-(u^{n+1}+u^n)_x/2
    Ify(kfy*ko2+1:kfy*ko2+ko2)=[Iy1l, Iy1l+1, Iy1l+2]+6*ko;
    Jfy(kfy*ko2+1:kfy*ko2+ko2)=[Jy1r, Jy1r, Jy1r]+ko;
    sfy(kfy*ko2+1:kfy*ko2+ko2)=-0.5*[s2, 3*s2, 3*s12xi]*h1e;  kfy=kfy+1;
    Ify(kfy*ko2+1:kfy*ko2+ko2)=[Iy1r, Iy1r+1, Iy1r+2]+6*ko;
    Jfy(kfy*ko2+1:kfy*ko2+ko2)=[Jy1r, Jy1r, Jy1r]+ko;
    sfy(kfy*ko2+1:kfy*ko2+ko2)=-0.5*[-s2, 3*s2, -3*s2xi]*h2e; kfy=kfy+1; 
    %%%%%%%%%%%%%%%%%%%% q5-(v^{n+1}+v^n)_x/2
    Ify(kfy*ko2+1:kfy*ko2+ko2)=[Iy1l, Iy1l+1, Iy1l+2]+8*ko;
    Jfy(kfy*ko2+1:kfy*ko2+ko2)=[Jy1r, Jy1r, Jy1r]+2*ko;
    sfy(kfy*ko2+1:kfy*ko2+ko2)=-0.5*[s2, 3*s2, 3*s12xi]*h1e;  kfy=kfy+1;
    Ify(kfy*ko2+1:kfy*ko2+ko2)=[Iy1r, Iy1r+1, Iy1r+2]+8*ko;
    Jfy(kfy*ko2+1:kfy*ko2+ko2)=[Jy1r, Jy1r, Jy1r]+2*ko;
    sfy(kfy*ko2+1:kfy*ko2+ko2)=-0.5*[-s2, 3*s2, -3*s2xi]*h2e; kfy=kfy+1;  
end

for k=1:fExn 
    id=fEx(k);   EK=Ex(id);   K1=Cells(EK.K1);    K2=Cells(EK.K2);   
    K1it=Cellsit(EK.K1);    K2it=Cellsit(EK.K2);  h1=K1.Wid;  h2=K2.Wid;  h=min(h1, h2);
    EK1f=K1.Indx_fine;  EK2f=K2.Indx_fine;  h1e=h/h1/h1; h2e=h/h2/h2;
    if K1.Lev==K2.Lev
        xib1=gs5;   UBn=K1.U*gpB;   UB=K1it.U*gpB;  
        xib2=gs5;   UTn=K2.U*gpT;   UT=K2it.U*gpT;    
    elseif K1.Lev==K2.Lev-1
        if K1.Center(1)>K2.Center(1)
            xib1=gps1;  UBn=K1.U*gcB1;   UB=K1it.U*gcB1; 
        else
            xib1=gps2;  UBn=K1.U*gcB2;   UB=K1it.U*gcB2; 
        end
        UTn=K2.U*gpT;    xib2=gs5;   UT=K2it.U*gpT; 
    elseif K1.Lev==K2.Lev+1
        if K2.Center(1)>K1.Center(1)
           UTn=K2.U*gcT1;   xib2=gps1;   UT=K2it.U*gcT1;
        else
           UTn=K2.U*gcT2;   xib2=gps2;   UT=K2it.U*gcT2;
        end
        xib1=gs5;   UBn=K1.U*gpB;  UB=K1it.U*gpB;
    else
        disp('Init: wrong conforming')
    end 
    rTn=UTn(1,:);  vTn=UTn(3,:);   rT=UT(1,:);  vT=UT(3,:);  zT=UT(4,:);
    uBn=UBn(2,:);  uTn=UTn(2,:);   uB=UB(2,:);  uT=UT(2,:);  zB=UB(4,:);
    nEK1 =(EK1f-1)*nk;   nEK2 =(EK2f-1)*nk;
    Iy1l=(nEK1+1)*ones(1,ko);   Iy1r=(nEK2+1)*ones(1,ko);   
    Jy1l=nEK1+1:nEK1+ko;    Jy1r=nEK2+1:nEK2+ko;  
%     %%%%%%%%%%%%%%%%%%%%%%%%r_t+(r*v)_y%%%%%%%%%%%%%%%%%%%
    s2e=(vT+vTn)/4;  
    s2   =[s2e*wt5, s2e.*xib2*wt5, -s2e*wt5]*0.5;
    s12xi=[s2e.*xib1*wt5, s2e.*xib2.*xib1*wt5, -s2e.*xib1*wt5]*0.5;
    s2xi =[s2e.*xib2*wt5, s2e.*xib2.*xib2*wt5, -s2e.*xib2*wt5]*0.5;
    Ify(kfy*ko2+1:kfy*ko2+ko2)=[Iy1l, Iy1l+1, Iy1l+2];
    Jfy(kfy*ko2+1:kfy*ko2+ko2)=[Jy1r, Jy1r, Jy1r];
    sfy(kfy*ko2+1:kfy*ko2+ko2)=dt*[s2, 3*s12xi, 3*s2]*h1e;  kfy=kfy+1;
    Ify(kfy*ko2+1:kfy*ko2+ko2)=[Iy1r, Iy1r+1, Iy1r+2];
    Jfy(kfy*ko2+1:kfy*ko2+ko2)=[Jy1r, Jy1r, Jy1r];
    sfy(kfy*ko2+1:kfy*ko2+ko2)=dt*[-s2, -3*s2xi, 3*s2]*h2e; kfy=kfy+1;  
    s2e=(rT+rTn)/4;  
    s2   =[s2e*wt5, s2e.*xib2*wt5, -s2e*wt5]*0.5;
    s12xi=[s2e.*xib1*wt5, s2e.*xib2.*xib1*wt5, -s2e.*xib1*wt5]*0.5;
    s2xi =[s2e.*xib2*wt5, s2e.*xib2.*xib2*wt5, -s2e.*xib2*wt5]*0.5;
    Ify(kfy*ko2+1:kfy*ko2+ko2)=[Iy1l, Iy1l+1, Iy1l+2];
    Jfy(kfy*ko2+1:kfy*ko2+ko2)=[Jy1r, Jy1r, Jy1r]+2*ko;
    sfy(kfy*ko2+1:kfy*ko2+ko2)=dt*[s2, 3*s12xi, 3*s2]*h1e;  kfy=kfy+1;
    Ify(kfy*ko2+1:kfy*ko2+ko2)=[Iy1r, Iy1r+1, Iy1r+2];
    Jfy(kfy*ko2+1:kfy*ko2+ko2)=[Jy1r, Jy1r, Jy1r]+2*ko;
    sfy(kfy*ko2+1:kfy*ko2+ko2)=dt*[-s2, -3*s2xi, 3*s2]*h2e; kfy=kfy+1;  
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     rZu=(rT+rTn)/8.*(vT+vTn).*(uT+uTn-uB-uBn);
%     F2f(EK2f,:)=F2f(EK2f,:)+[rZu*wt5, rZu.*xib2*wt5,-rZu*wt5]*0.5*h; 
    s2e=1/8.*(vT+vTn).*(uT+uTn-uB-uBn);  
    s2   =[s2e*wt5, s2e.*xib2*wt5, -s2e*wt5]*0.5;
    s2xi =[s2e.*xib2*wt5, s2e.*xib2.*xib2*wt5, -s2e.*xib2*wt5]*0.5;
    Ify(kfy*ko2+1:kfy*ko2+ko2)=[Iy1r, Iy1r+1, Iy1r+2]+ko;
    Jfy(kfy*ko2+1:kfy*ko2+ko2)=[Jy1r, Jy1r, Jy1r];
    sfy(kfy*ko2+1:kfy*ko2+ko2)=dt*[s2, s2xi, -s2]*h2e; kfy=kfy+1; 
    
    s2e=1/8.*(rT+rTn).*(uT+uTn-uB-uBn);  
    s2   =[s2e*wt5, s2e.*xib2*wt5, -s2e*wt5]*0.5;
    s2xi =[s2e.*xib2*wt5, s2e.*xib2.*xib2*wt5, -s2e.*xib2*wt5]*0.5;
    Ify(kfy*ko2+1:kfy*ko2+ko2)=[Iy1r, Iy1r+1, Iy1r+2]+ko;
    Jfy(kfy*ko2+1:kfy*ko2+ko2)=[Jy1r, Jy1r, Jy1r]+2*ko;
    sfy(kfy*ko2+1:kfy*ko2+ko2)=dt*[s2, s2xi, -s2]*h2e; kfy=kfy+1; 
    
    s2e=1/8*(rT+rTn).*(vT+vTn);  
    s2   =[s2e*wt5, s2e.*xib2*wt5, -s2e*wt5]*0.5;
    s2xi =[s2e.*xib2*wt5, s2e.*xib2.*xib2*wt5, -s2e.*xib2*wt5]*0.5;
    Ify(kfy*ko2+1:kfy*ko2+ko2)=[Iy1r, Iy1r+1, Iy1r+2]+ko;
    Jfy(kfy*ko2+1:kfy*ko2+ko2)=[Jy1r, Jy1r, Jy1r]+ko;
    sfy(kfy*ko2+1:kfy*ko2+ko2)=dt*[s2, s2xi, -s2]*h2e; kfy=kfy+1; 
    
    s1e=-1/8*(rT+rTn).*(vT+vTn); 
    s1   =[s1e*wt5, s1e.*xib1*wt5, s1e*wt5]*0.5;  
    s21xi=[s1e.*xib2*wt5, s1e.*xib1.*xib2*wt5, s1e.*xib2*wt5]*0.5;
    Ify(kfy*ko2+1:kfy*ko2+ko2)=[Iy1r, Iy1r+1, Iy1r+2]+ko;
    Jfy(kfy*ko2+1:kfy*ko2+ko2)=[Jy1l, Jy1l, Jy1l]+ko;
    sfy(kfy*ko2+1:kfy*ko2+ko2)=dt*[s1, s21xi, -s1]*h2e; kfy=kfy+1;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     rZu=(rT+rTn)/2.*((zT-uT.^2/4-uTn.^2/4)-(zB-uB.^2/4-uBn.^2/4));
%     F3f(EK2f,:)=F3f(EK2f,:)+[rZu*wt5, rZu.*xib2*wt5,-rZu*wt5]*0.5*h; 
    s2e=1/2.*((zT-uT.^2/4-uTn.^2/4)-(zB-uB.^2/4-uBn.^2/4));
    s2   =[s2e*wt5, s2e.*xib2*wt5, -s2e*wt5]*0.5;
    s2xi =[s2e.*xib2*wt5, s2e.*xib2.*xib2*wt5, -s2e.*xib2*wt5]*0.5;
    Ify(kfy*ko2+1:kfy*ko2+ko2)=[Iy1r, Iy1r+1, Iy1r+2]+2*ko;
    Jfy(kfy*ko2+1:kfy*ko2+ko2)=[Jy1r, Jy1r, Jy1r];
    sfy(kfy*ko2+1:kfy*ko2+ko2)=dt*[s2, s2xi, -s2]*h2e; kfy=kfy+1; 
    s2e=(rT+rTn)/2;
    s2   =[s2e*wt5, s2e.*xib2*wt5, -s2e*wt5]*0.5;
    s2xi =[s2e.*xib2*wt5, s2e.*xib2.*xib2*wt5, -s2e.*xib2*wt5]*0.5;
    Ify(kfy*ko2+1:kfy*ko2+ko2)=[Iy1r, Iy1r+1, Iy1r+2]+2*ko;
    Jfy(kfy*ko2+1:kfy*ko2+ko2)=[Jy1r, Jy1r, Jy1r]+3*ko;
    sfy(kfy*ko2+1:kfy*ko2+ko2)=dt*[s2, s2xi, -s2]*h2e; kfy=kfy+1; 
    s1e=-s2e;
    s1   =[s1e*wt5, s1e.*xib1*wt5, s1e*wt5]*0.5;  
    s21xi=[s1e.*xib2*wt5, s1e.*xib1.*xib2*wt5, s1e.*xib2*wt5]*0.5;
    Ify(kfy*ko2+1:kfy*ko2+ko2)=[Iy1r, Iy1r+1, Iy1r+2]+2*ko;
    Jfy(kfy*ko2+1:kfy*ko2+ko2)=[Jy1l, Jy1l, Jy1l]+3*ko;
    sfy(kfy*ko2+1:kfy*ko2+ko2)=dt*[s1, s21xi,-s1]*h2e; kfy=kfy+1;
    
    s2e=(rT+rTn)/2.*(-uT/2);
    s2   =[s2e*wt5, s2e.*xib2*wt5, -s2e*wt5]*0.5;
    s2xi =[s2e.*xib2*wt5, s2e.*xib2.*xib2*wt5, -s2e.*xib2*wt5]*0.5;
    Ify(kfy*ko2+1:kfy*ko2+ko2)=[Iy1r, Iy1r+1, Iy1r+2]+2*ko;
    Jfy(kfy*ko2+1:kfy*ko2+ko2)=[Jy1r, Jy1r, Jy1r]+ko;
    sfy(kfy*ko2+1:kfy*ko2+ko2)=dt*[s2, s2xi, -s2]*h2e; kfy=kfy+1; 
    s1e=(rT+rTn)/2.*(uB/2);  
    s1   =[s1e*wt5, s1e.*xib1*wt5, s1e*wt5]*0.5;  
    s21xi=[s1e.*xib2*wt5, s1e.*xib1.*xib2*wt5, s1e.*xib2*wt5]*0.5;
    Ify(kfy*ko2+1:kfy*ko2+ko2)=[Iy1r, Iy1r+1, Iy1r+2]+2*ko;
    Jfy(kfy*ko2+1:kfy*ko2+ko2)=[Jy1l, Jy1l, Jy1l]+ko;
    sfy(kfy*ko2+1:kfy*ko2+ko2)=dt*[s1, s21xi,-s1]*h2e; kfy=kfy+1;
    
%     %%%%%%%%%%%%%%%%%%%%%%% Z+lambda*q2_y
    s1e=oneb;  
    s1   =[s1e*wt5, s1e.*xib1*wt5, s1e*wt5]*0.5;  
    s1xi =[s1e.*xib1*wt5, s1e.*xib1.*xib1*wt5, s1e.*xib1*wt5]*0.5;
    s21xi=[s1e.*xib2*wt5, s1e.*xib1.*xib2*wt5, s1e.*xib2*wt5]*0.5;
%     %     f2y=-mu*(Ukxi(8,:)+Ukxi(9,:));
    Ify(kfy*ko2+1:kfy*ko2+ko2)=[Iy1l, Iy1l+1, Iy1l+2]+ko;
    Jfy(kfy*ko2+1:kfy*ko2+ko2)=[Jy1l, Jy1l, Jy1l]+7*ko;
    sfy(kfy*ko2+1:kfy*ko2+ko2)=-mu*dt*[s1, s1xi, s1]*h1e;   kfy=kfy+1;
    Ify(kfy*ko2+1:kfy*ko2+ko2)=[Iy1r, Iy1r+1, Iy1r+2]+ko;
    Jfy(kfy*ko2+1:kfy*ko2+ko2)=[Jy1l, Jy1l, Jy1l]+7*ko;
    sfy(kfy*ko2+1:kfy*ko2+ko2)=-mu*dt*[-s1,-s21xi, s1]*h2e; kfy=kfy+1;
    
    Ify(kfy*ko2+1:kfy*ko2+ko2)=[Iy1l, Iy1l+1, Iy1l+2]+ko;
    Jfy(kfy*ko2+1:kfy*ko2+ko2)=[Jy1l, Jy1l, Jy1l]+8*ko;
    sfy(kfy*ko2+1:kfy*ko2+ko2)=-mu*dt*[s1, s1xi, s1]*h1e;   kfy=kfy+1;
    Ify(kfy*ko2+1:kfy*ko2+ko2)=[Iy1r, Iy1r+1, Iy1r+2]+ko;
    Jfy(kfy*ko2+1:kfy*ko2+ko2)=[Jy1l, Jy1l, Jy1l]+8*ko;
    sfy(kfy*ko2+1:kfy*ko2+ko2)=-mu*dt*[-s1,-s21xi, s1]*h2e; kfy=kfy+1;
%      f3y=-mu*(4/3*Ukxi(10,:)-2/3*Ukxi(7,:)); 
        Ify(kfy*ko2+1:kfy*ko2+ko2)=[Iy1l, Iy1l+1, Iy1l+2]+2*ko;
    Jfy(kfy*ko2+1:kfy*ko2+ko2)=[Jy1l, Jy1l, Jy1l]+9*ko;
    sfy(kfy*ko2+1:kfy*ko2+ko2)=-4/3*mu*dt*[s1, s1xi, s1]*h1e;   kfy=kfy+1;
    Ify(kfy*ko2+1:kfy*ko2+ko2)=[Iy1r, Iy1r+1, Iy1r+2]+2*ko;
    Jfy(kfy*ko2+1:kfy*ko2+ko2)=[Jy1l, Jy1l, Jy1l]+9*ko;
    sfy(kfy*ko2+1:kfy*ko2+ko2)=-4/3*mu*dt*[-s1,-s21xi, s1]*h2e; kfy=kfy+1;
    
    Ify(kfy*ko2+1:kfy*ko2+ko2)=[Iy1l, Iy1l+1, Iy1l+2]+2*ko;
    Jfy(kfy*ko2+1:kfy*ko2+ko2)=[Jy1l, Jy1l, Jy1l]+6*ko;
    sfy(kfy*ko2+1:kfy*ko2+ko2)=2/3*mu*dt*[s1, s1xi, s1]*h1e;   kfy=kfy+1;
    Ify(kfy*ko2+1:kfy*ko2+ko2)=[Iy1r, Iy1r+1, Iy1r+2]+2*ko;
    Jfy(kfy*ko2+1:kfy*ko2+ko2)=[Jy1l, Jy1l, Jy1l]+6*ko;
    sfy(kfy*ko2+1:kfy*ko2+ko2)=2/3*mu*dt*[-s1,-s21xi, s1]*h2e; kfy=kfy+1;
%     %%%%%%%%%%%%%%%%%%%%% Z
    Ify(kfy*ko2+1:kfy*ko2+ko2)=[Iy1l, Iy1l+1, Iy1l+2]+3*ko;
    Jfy(kfy*ko2+1:kfy*ko2+ko2)=[Jy1l, Jy1l, Jy1l]+5*ko;
    sfy(kfy*ko2+1:kfy*ko2+ko2)=lambda*[s1, 3.0*s1xi, 3.0*s1]*h1e;   kfy=kfy+1;
    Ify(kfy*ko2+1:kfy*ko2+ko2)=[Iy1r, Iy1r+1, Iy1r+2]+3*ko;
    Jfy(kfy*ko2+1:kfy*ko2+ko2)=[Jy1l, Jy1l, Jy1l]+5*ko;
    sfy(kfy*ko2+1:kfy*ko2+ko2)=lambda*[-s1,-3.0*s21xi, 3.0*s1]*h2e; kfy=kfy+1;
    %%%%%%%%%%%%%%%%%%%%%% q2=r_y
    s2e=oneb;
    s2   =[s2e*wt5, s2e.*xib2*wt5, -s2e*wt5]*0.5;
    s12xi=[s2e.*xib1*wt5, s2e.*xib2.*xib1*wt5, -s2e.*xib1*wt5]*0.5;
    s2xi =[s2e.*xib2*wt5, s2e.*xib2.*xib2*wt5, -s2e.*xib2*wt5]*0.5;
    Ify(kfy*ko2+1:kfy*ko2+ko2)=[Iy1l, Iy1l+1, Iy1l+2]+5*ko;
    Jfy(kfy*ko2+1:kfy*ko2+ko2)=[Jy1r, Jy1r, Jy1r];
    sfy(kfy*ko2+1:kfy*ko2+ko2)=-0.5*[s2, 3*s12xi, 3*s2]*h1e;  kfy=kfy+1;
    Ify(kfy*ko2+1:kfy*ko2+ko2)=[Iy1r, Iy1r+1, Iy1r+2]+5*ko;
    Jfy(kfy*ko2+1:kfy*ko2+ko2)=[Jy1r, Jy1r, Jy1r];
    sfy(kfy*ko2+1:kfy*ko2+ko2)=-0.5*[-s2, -3*s2xi, 3*s2]*h2e; kfy=kfy+1;  
    %%%%%%%%%%%%%%%%%%%% q4-(u^{n+1}+u^n)_y/2
    Ify(kfy*ko2+1:kfy*ko2+ko2)=[Iy1l, Iy1l+1, Iy1l+2]+7*ko;
    Jfy(kfy*ko2+1:kfy*ko2+ko2)=[Jy1r, Jy1r, Jy1r]+ko;
    sfy(kfy*ko2+1:kfy*ko2+ko2)=-0.5*[s2, 3*s12xi, 3*s2]*h1e;  kfy=kfy+1;
    Ify(kfy*ko2+1:kfy*ko2+ko2)=[Iy1r, Iy1r+1, Iy1r+2]+7*ko;
    Jfy(kfy*ko2+1:kfy*ko2+ko2)=[Jy1r, Jy1r, Jy1r]+ko;
    sfy(kfy*ko2+1:kfy*ko2+ko2)=-0.5*[-s2, -3*s2xi, 3*s2]*h2e; kfy=kfy+1;  
    %%%%%%%%%%%%%%%%%%%% q6-(v^{n+1}+v^n)_y/2
    Ify(kfy*ko2+1:kfy*ko2+ko2)=[Iy1l, Iy1l+1, Iy1l+2]+9*ko;
    Jfy(kfy*ko2+1:kfy*ko2+ko2)=[Jy1r, Jy1r, Jy1r]+2*ko;
    sfy(kfy*ko2+1:kfy*ko2+ko2)=-0.5*[s2, 3*s12xi, 3*s2]*h1e;  kfy=kfy+1;
    Ify(kfy*ko2+1:kfy*ko2+ko2)=[Iy1r, Iy1r+1, Iy1r+2]+9*ko;
    Jfy(kfy*ko2+1:kfy*ko2+ko2)=[Jy1r, Jy1r, Jy1r]+2*ko;
    sfy(kfy*ko2+1:kfy*ko2+ko2)=-0.5*[-s2, -3*s2xi, 3*s2]*h2e; kfy=kfy+1;  
end

onek=ones(1,gk*gk);  dek=[onek*wt, onek.*xxi*wt, onek.*yyi*wt];
Mrut_r=zeros(ko,ko);    Mrut_u=zeros(ko,ko);    Mrvt_r=zeros(ko,ko); 
MZ_u=zeros(ko,ko);  MZ_v=MZ_u;  MZ_r=MZ_u;
for k=1:fCellsn
    id=fCells(k);   Ckn=Cells(id);  h=Ckn.Wid;  Ukn=Ckn.U;  Uknxi=Ukn*gpM;
    rxin=Uknxi(1,:);  uxin=Uknxi(2,:);   vxin=Uknxi(3,:);  
    Uk=Cellsit(id).U;    Ukxi=Uk*gpM; 
    rxi=Ukxi(1,:);  uxi=Ukxi(2,:);   vxi=Ukxi(3,:);
%     rxi=Ukxi(1,:);   ritxi=Ukitxi(1,:);  
    IJy1=(k-1)*nk+1:(k-1)*nk+ko; 
    rkc=(rxi+rxin)/2;  vkc=(vxi+vxin)/2;    ukc=(uxi+uxin)/2;
    Mruc_u=zeros(ko,ko);  Mruc_v=zeros(ko,ko);    MruZ_r=zeros(ko,ko);    MruZ_Z=zeros(ko,ko);
    Mrvc_u=zeros(ko,ko);  Mrvc_v=zeros(ko,ko);    MrvZ_r=zeros(ko,ko);    MrvZ_Z=zeros(ko,ko);
    if norm(Uk(1,:)-Ukn(1,:), inf)<myeps
        Z_rk=Wrr((rxi+rxin)/2)/2;
    else
        Z_rk=(Wr(rxi).*(rxi-rxin)-(W(rxi)-W(rxin)))./(rxi-rxin).^2;
    end
    for j=1:ko
        for i=1:ko
            Mrut_r(i,j)=1/4.*(uxi-uxin)/2.*gpM(i,:).*gpM(j,:)*wt;    
            Mrut_u(i,j)=1/4.*(rxi+rxin)/2.*gpM(i,:).*gpM(j,:)*wt;     
            Mrvt_r(i,j)=1/4.*(vxi-vxin)/2.*gpM(i,:).*gpM(j,:)*wt;    
            MZ_u(i,j)=-1/4*aav(i)*1/2*uxi.*gpM(j,:).*gpM(i,:)*wt;
            MZ_v(i,j)=-1/4*aav(i)*1/2*vxi.*gpM(j,:).*gpM(i,:)*wt;
            MZ_r(i,j)=-1/4*aav(i)*(Z_rk).*gpM(j,:).*gpM(i,:)*wt;
            %%%%%%%%%%%%%%%%%%%%%%%
            ruc_v=rkc.*gpM(i,:).*(-gpM(j,:).*(Uk(3,:)*tfx)/2-vxi.*tfx(j,:)/2)*wt+...
                  rkc.*gpM(j,:)/2.*gpM(i,:).*(Uk(2,:)*tfy+Ukn(2,:)*tfy)/2*wt;
            Mruc_v(i,j)=dt/h/2*ruc_v;     
            rZu_r=1/2*gpM(j,:).*gpM(i,:).*(Uk(4,:)*tfx-vxi.*(Uk(3,:)*tfx)/2-vxin.*(Ukn(3,:)*tfx)/2)*wt...
                +1/2*gpM(j,:).*vkc.*gpM(i,:).*(Uk(2,:)*tfy+Ukn(2,:)*tfy)/2*wt;  
            MruZ_r(i,j)=dt/h/2*rZu_r;
            rZ_Zk=rkc.*gpM(i,:).*(tfx(j,:))*wt;    MruZ_Z(i,j)=dt/h/2*rZ_Zk;
            ruc_u=rkc.*vkc.*gpM(i,:).*tfy(j,:)/2*wt; Mruc_u(i,j)=dt/h/2*ruc_u;
            %%%%%%%%%%%%%%%%%%%%%%
            rvc_u=rkc.*gpM(i,:).*(-gpM(j,:).*(Uk(2,:)*tfy)/2-uxi.*tfy(j,:)/2)*wt+...
                rkc.*gpM(j,:)/2.*gpM(i,:).*(Uk(3,:)*tfx+Ukn(3,:)*tfx)/2*wt;
            Mrvc_u(i,j)=dt/h/2*rvc_u;
            rvc_v=rkc.*ukc.*gpM(i,:).*tfx(j,:)/2*wt;      Mrvc_v(i,j)=dt/h/2*rvc_v;     
    
            rZv_r=1/2*gpM(j,:).*gpM(i,:).*(Uk(4,:)*tfy-uxi.*(Uk(2,:)*tfy)/2-uxin.*(Ukn(2,:)*tfy)/2)*wt+...
                 1/2*gpM(j,:).*ukc.*gpM(i,:).*(Uk(3,:)*tfx+Ukn(3,:)*tfx)/2*wt;  
            MrvZ_r(i,j)=dt/h/2*rZv_r;
            rZ_Zk=rkc.*gpM(i,:).*(tfy(j,:))*wt;   MrvZ_Z(i,j)=dt/h/2*rZ_Zk;
        end
    end
    Ify(kfy*ko2+1:kfy*ko2+ko2)=[IJy1, IJy1, IJy1]+ko;
    Jfy(kfy*ko2+1:kfy*ko2+ko2)=[IJy1; IJy1; IJy1]+ko;
    sfy(kfy*ko2+1:kfy*ko2+ko2)=Mruc_u;   kfy=kfy+1;
    Ify(kfy*ko2+1:kfy*ko2+ko2)=[IJy1, IJy1, IJy1]+ko;
    Jfy(kfy*ko2+1:kfy*ko2+ko2)=[IJy1; IJy1; IJy1]+2*ko;
    sfy(kfy*ko2+1:kfy*ko2+ko2)=Mruc_v;   kfy=kfy+1;
    Ify(kfy*ko2+1:kfy*ko2+ko2)=[IJy1, IJy1, IJy1]+ko;
    Jfy(kfy*ko2+1:kfy*ko2+ko2)=[IJy1; IJy1; IJy1]+3*ko;
    sfy(kfy*ko2+1:kfy*ko2+ko2)=MruZ_Z;   kfy=kfy+1;
    Ify(kfy*ko2+1:kfy*ko2+ko2)=[IJy1, IJy1, IJy1];
    Jfy(kfy*ko2+1:kfy*ko2+ko2)=[IJy1; IJy1; IJy1];
    sfy(kfy*ko2+1:kfy*ko2+ko2)=eye(ko); kfy=kfy+1; 
    Ify(kfy*ko2+1:kfy*ko2+ko2)=[IJy1, IJy1, IJy1]+ko;
    Jfy(kfy*ko2+1:kfy*ko2+ko2)=[IJy1; IJy1; IJy1]+ko;
    sfy(kfy*ko2+1:kfy*ko2+ko2)=Mrut_u;   kfy=kfy+1;
    Ify(kfy*ko2+1:kfy*ko2+ko2)=[IJy1, IJy1, IJy1]+ko;
    Jfy(kfy*ko2+1:kfy*ko2+ko2)=[IJy1; IJy1; IJy1];
    sfy(kfy*ko2+1:kfy*ko2+ko2)=Mrut_r+MruZ_r;   kfy=kfy+1;%%
   
    Ify(kfy*ko2+1:kfy*ko2+ko2)=[IJy1, IJy1, IJy1]+2*ko;
    Jfy(kfy*ko2+1:kfy*ko2+ko2)=[IJy1; IJy1; IJy1]+ko;
    sfy(kfy*ko2+1:kfy*ko2+ko2)=Mrvc_u;   kfy=kfy+1;
    Ify(kfy*ko2+1:kfy*ko2+ko2)=[IJy1, IJy1, IJy1]+2*ko;
    Jfy(kfy*ko2+1:kfy*ko2+ko2)=[IJy1; IJy1; IJy1]+2*ko;
    sfy(kfy*ko2+1:kfy*ko2+ko2)=Mrvc_v;   kfy=kfy+1;
    Ify(kfy*ko2+1:kfy*ko2+ko2)=[IJy1, IJy1, IJy1]+2*ko;
    Jfy(kfy*ko2+1:kfy*ko2+ko2)=[IJy1; IJy1; IJy1]+3*ko;
    sfy(kfy*ko2+1:kfy*ko2+ko2)=MrvZ_Z;   kfy=kfy+1;
    Ify(kfy*ko2+1:kfy*ko2+ko2)=[IJy1, IJy1, IJy1]+2*ko;
    Jfy(kfy*ko2+1:kfy*ko2+ko2)=[IJy1; IJy1; IJy1];
    sfy(kfy*ko2+1:kfy*ko2+ko2)=Mrvt_r+MrvZ_r;   kfy=kfy+1;
    Ify(kfy*ko2+1:kfy*ko2+ko2)=[IJy1, IJy1, IJy1]+2*ko;
    Jfy(kfy*ko2+1:kfy*ko2+ko2)=[IJy1; IJy1; IJy1]+2*ko;
    sfy(kfy*ko2+1:kfy*ko2+ko2)=Mrut_u;   kfy=kfy+1;
    Ify(kfy*ko2+1:kfy*ko2+ko2)=[IJy1, IJy1, IJy1]+3*ko;
    Jfy(kfy*ko2+1:kfy*ko2+ko2)=[IJy1; IJy1; IJy1]+3*ko;
    sfy(kfy*ko2+1:kfy*ko2+ko2)=eye(ko);  kfy=kfy+1; 
    Ify(kfy*ko2+1:kfy*ko2+ko2)=[IJy1, IJy1, IJy1]+3*ko;
    Jfy(kfy*ko2+1:kfy*ko2+ko2)=[IJy1; IJy1; IJy1];
    sfy(kfy*ko2+1:kfy*ko2+ko2)=MZ_r; kfy=kfy+1; 
    Ify(kfy*ko2+1:kfy*ko2+ko2)=[IJy1, IJy1, IJy1]+3*ko;
    Jfy(kfy*ko2+1:kfy*ko2+ko2)=[IJy1; IJy1; IJy1]+ko;
    sfy(kfy*ko2+1:kfy*ko2+ko2)=MZ_u;    kfy=kfy+1; 
    Ify(kfy*ko2+1:kfy*ko2+ko2)=[IJy1, IJy1, IJy1]+3*ko;
    Jfy(kfy*ko2+1:kfy*ko2+ko2)=[IJy1; IJy1; IJy1]+2*ko;
    sfy(kfy*ko2+1:kfy*ko2+ko2)=MZ_v; kfy=kfy+1; 
    
    Ify(kfy*ko2+1:kfy*ko2+ko2)=[IJy1, IJy1, IJy1]+4*ko;
    Jfy(kfy*ko2+1:kfy*ko2+ko2)=[IJy1; IJy1; IJy1]+4*ko;
    sfy(kfy*ko2+1:kfy*ko2+ko2)=eye(ko); kfy=kfy+1; 
    Ify(kfy*ko2+1:kfy*ko2+ko2)=[IJy1, IJy1, IJy1]+5*ko;
    Jfy(kfy*ko2+1:kfy*ko2+ko2)=[IJy1; IJy1; IJy1]+5*ko;
    sfy(kfy*ko2+1:kfy*ko2+ko2)=eye(ko); kfy=kfy+1; 
    Ify(kfy*ko2+1:kfy*ko2+ko2)=[IJy1, IJy1, IJy1]+6*ko;
    Jfy(kfy*ko2+1:kfy*ko2+ko2)=[IJy1; IJy1; IJy1]+6*ko;
    sfy(kfy*ko2+1:kfy*ko2+ko2)=eye(ko); kfy=kfy+1; 
    Ify(kfy*ko2+1:kfy*ko2+ko2)=[IJy1, IJy1, IJy1]+7*ko;
    Jfy(kfy*ko2+1:kfy*ko2+ko2)=[IJy1; IJy1; IJy1]+7*ko;
    sfy(kfy*ko2+1:kfy*ko2+ko2)=eye(ko); kfy=kfy+1; 
    Ify(kfy*ko2+1:kfy*ko2+ko2)=[IJy1, IJy1, IJy1]+8*ko;
    Jfy(kfy*ko2+1:kfy*ko2+ko2)=[IJy1; IJy1; IJy1]+8*ko;
    sfy(kfy*ko2+1:kfy*ko2+ko2)=eye(ko); kfy=kfy+1; 
    Ify(kfy*ko2+1:kfy*ko2+ko2)=[IJy1, IJy1, IJy1]+9*ko;
    Jfy(kfy*ko2+1:kfy*ko2+ko2)=[IJy1; IJy1; IJy1]+9*ko;
    sfy(kfy*ko2+1:kfy*ko2+ko2)=eye(ko); kfy=kfy+1; 
%     %%%%%%%%%%%%%%% r_t+(r*u)_x+(r*v)_y
    uck=(uxi+uxin)/4;  
    Ifye(kfye*ko+1:kfye*ko+ko)=IJy1(1)+1;   Jfye(kfye*ko+1:kfye*ko+ko)=IJy1;
    sfye(kfye*ko+1:kfye*ko+ko)=-dt*[uck*wt, uck.*xxi*wt, uck.*yyi*wt]*1.5/h;  kfye=kfye+1;
    vck=(vxi+vxin)/4;  
    Ifye(kfye*ko+1:kfye*ko+ko)=IJy1(1)+2;   Jfye(kfye*ko+1:kfye*ko+ko)=IJy1;
    sfye(kfye*ko+1:kfye*ko+ko)=-dt*[vck*wt, vck.*xxi*wt, vck.*yyi*wt]*1.5/h;  kfye=kfye+1; 
    rck=(rxi+rxin)/4;  
    Ifye(kfye*ko+1:kfye*ko+ko)=IJy1(1)+1;   Jfye(kfye*ko+1:kfye*ko+ko)=IJy1+ko;
    sfye(kfye*ko+1:kfye*ko+ko)=-dt*[rck*wt, rck.*xxi*wt, rck.*yyi*wt]*1.5/h;  kfye=kfye+1; 
    Ifye(kfye*ko+1:kfye*ko+ko)=IJy1(1)+2;   Jfye(kfye*ko+1:kfye*ko+ko)=IJy1+2*ko;
    sfye(kfye*ko+1:kfye*ko+ko)=-dt*[rck*wt, rck.*xxi*wt, rck.*yyi*wt]*1.5/h;  kfye=kfye+1;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %%%%%%%%%%%%%%%%%%%%  
%  f2x=-mu*(4/3*Ukxi(7,:)-2/3*Ukxi(10,:));  f2y=-mu*(Ukxi(8,:)+Ukxi(9,:));
    Ifye(kfye*ko+1:kfye*ko+ko)=IJy1(1)+ko+1;     Jfye(kfye*ko+1:kfye*ko+ko)=IJy1+6*ko;
    sfye(kfye*ko+1:kfye*ko+ko)=4/3*mu*dt*dek*0.5/h;  kfye=kfye+1;
    Ifye(kfye*ko+1:kfye*ko+ko)=IJy1(1)+ko+1;     Jfye(kfye*ko+1:kfye*ko+ko)=IJy1+9*ko;
    sfye(kfye*ko+1:kfye*ko+ko)=-2/3*mu*dt*dek*0.5/h;  kfye=kfye+1;
    Ifye(kfye*ko+1:kfye*ko+ko)=IJy1(1)+ko+2;     Jfye(kfye*ko+1:kfye*ko+ko)=IJy1+7*ko;
    sfye(kfye*ko+1:kfye*ko+ko)=mu*dt*dek*0.5/h;  kfye=kfye+1;
    Ifye(kfye*ko+1:kfye*ko+ko)=IJy1(1)+ko+2;     Jfye(kfye*ko+1:kfye*ko+ko)=IJy1+8*ko;
    sfye(kfye*ko+1:kfye*ko+ko)=mu*dt*dek*0.5/h;  kfye=kfye+1;
    %  f3y=-mu*(4/3*Ukxi(10,:)-2/3*Ukxi(7,:));  f3x=-mu*(Ukxi(8,:)+Ukxi(9,:));
    Ifye(kfye*ko+1:kfye*ko+ko)=IJy1(1)+2*ko+1;     Jfye(kfye*ko+1:kfye*ko+ko)=IJy1+7*ko;
    sfye(kfye*ko+1:kfye*ko+ko)=mu*dt*dek*0.5/h;  kfye=kfye+1;
    Ifye(kfye*ko+1:kfye*ko+ko)=IJy1(1)+2*ko+1;     Jfye(kfye*ko+1:kfye*ko+ko)=IJy1+8*ko;
    sfye(kfye*ko+1:kfye*ko+ko)=mu*dt*dek*0.5/h;  kfye=kfye+1;
    Ifye(kfye*ko+1:kfye*ko+ko)=IJy1(1)+2*ko+2;     Jfye(kfye*ko+1:kfye*ko+ko)=IJy1+6*ko;
    sfye(kfye*ko+1:kfye*ko+ko)=-2/3*mu*dt*dek*0.5/h;  kfye=kfye+1;
    Ifye(kfye*ko+1:kfye*ko+ko)=IJy1(1)+2*ko+2;     Jfye(kfye*ko+1:kfye*ko+ko)=IJy1+9*ko;
    sfye(kfye*ko+1:kfye*ko+ko)=4/3*mu*dt*dek*0.5/h;  kfye=kfye+1;
    %%%%%%%%%%%%%%%% Z+lambda*q1_x
    Ifye(kfye*ko+1:kfye*ko+ko)=IJy1(1)+3*ko+1;     Jfye(kfye*ko+1:kfye*ko+ko)=IJy1+4*ko;
    sfye(kfye*ko+1:kfye*ko+ko)=-lambda*dek*1.5/h;  kfye=kfye+1;
    %%%%%%%%%%%%%%%% Z+lambda*q2_y
    Ifye(kfye*ko+1:kfye*ko+ko)=IJy1(1)+3*ko+2;     Jfye(kfye*ko+1:kfye*ko+ko)=IJy1+5*ko;
    sfye(kfye*ko+1:kfye*ko+ko)=-lambda*dek*1.5/h;  kfye=kfye+1; 
    %%%%%%%%%%%%%%%% q1-r_x
    Ifye(kfye*ko+1:kfye*ko+ko)=IJy1(1)+4*ko+1;     Jfye(kfye*ko+1:kfye*ko+ko)=IJy1;
    sfye(kfye*ko+1:kfye*ko+ko)=0.5*dek*1.5/h;  kfye=kfye+1;
    %%%%%%%%%%%%%%%% q2-r_y
    Ifye(kfye*ko+1:kfye*ko+ko)=IJy1(1)+5*ko+2;     Jfye(kfye*ko+1:kfye*ko+ko)=IJy1;
    sfye(kfye*ko+1:kfye*ko+ko)=0.5*dek*1.5/h;  kfye=kfye+1; 
    %%%%%%%%%%%%%%%% q3-u_x
    Ifye(kfye*ko+1:kfye*ko+ko)=IJy1(1)+6*ko+1;     Jfye(kfye*ko+1:kfye*ko+ko)=IJy1+ko;
    sfye(kfye*ko+1:kfye*ko+ko)=0.5*dek*1.5/h;  kfye=kfye+1;
    %%%%%%%%%%%%%%%% q4-u_y
    Ifye(kfye*ko+1:kfye*ko+ko)=IJy1(1)+7*ko+2;     Jfye(kfye*ko+1:kfye*ko+ko)=IJy1+ko;
    sfye(kfye*ko+1:kfye*ko+ko)=0.5*dek*1.5/h;  kfye=kfye+1; 
    %%%%%%%%%%%%%%%% q5-v_x
    Ifye(kfye*ko+1:kfye*ko+ko)=IJy1(1)+8*ko+1;     Jfye(kfye*ko+1:kfye*ko+ko)=IJy1+2*ko;
    sfye(kfye*ko+1:kfye*ko+ko)=0.5*dek*1.5/h;  kfye=kfye+1;
    %%%%%%%%%%%%%%%% q6-v_y
    Ifye(kfye*ko+1:kfye*ko+ko)=IJy1(1)+9*ko+2;     Jfye(kfye*ko+1:kfye*ko+ko)=IJy1+2*ko;
    sfye(kfye*ko+1:kfye*ko+ko)=0.5*dek*1.5/h;  kfye=kfye+1; 
end
Mat=sparse(Ify(1:kfy*ko2), Jfy(1:kfy*ko2), sfy(1:kfy*ko2), nks, nks)+...
    sparse(Ifye(1:kfye*ko), Jfye(1:kfye*ko), sfye(1:kfye*ko), nks, nks); 
end