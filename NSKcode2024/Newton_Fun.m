function [Fu]=Newton_Fun(Cellsit, Mesh, Mesh_fine, t, dt)
global gs5 wt5 gps1 gps2 wt gsx gsy  mu lambda  Wr W tfx tfy myeps
global gpL gpR gpB gpT  gcL1 gcR1 gcB1 gcT1  gcL2 gcR2 gcB2 gcT2 gpM
ko=3;   nk=10*ko;  
%%%%%%%%%%%%%%%%%%%%% Mesh_fine u1,u2,p
fCells=Mesh_fine.fCells;    fCellsn=Mesh_fine.fCellsn;
fEy=Mesh_fine.fEy;   fEyn=Mesh_fine.fEyn;   
fEx=Mesh_fine.fEx;   fExn=Mesh_fine.fExn;             
Cells=Mesh.Cells;   Ey=Mesh.Ey;   Ex=Mesh.Ex;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
F1f=zeros(fCellsn,ko);  F2f=F1f;    F3f=F1f;  F4f=F1f;  F5f=F1f;   
F6f=F1f;  F7f=F1f;  F8f=F1f;    F9f=F1f;  F10f=F1f;
for k=1:fEyn 
    id=fEy(k);   EK=Ey(id);   K1=Cells(EK.K1);    K2=Cells(EK.K2);  
    K1it=Cellsit(EK.K1);    K2it=Cellsit(EK.K2);    
    EK1f=K1.Indx_fine;  EK2f=K2.Indx_fine;  h=min(K1.Wid, K2.Wid); 
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
    f1h=1/4*(rR.*uR+rR.*uRn+rRn.*uR+rRn.*uRn); 
    F1f(EK1f,:)=F1f(EK1f,:)+[f1h*wt5, f1h*wt5*3, f1h.*xib1*wt5*3]*0.5*h;  
    F1f(EK2f,:)=F1f(EK2f,:)-[f1h*wt5,-f1h*wt5*3, f1h.*xib2*wt5*3]*0.5*h; 
    qsh=-mu*(4/3*UL(7,:)-2/3*UL(10,:));
    F2f(EK1f,:)=F2f(EK1f,:)+[qsh*wt5, qsh*wt5, qsh.*xib1*wt5]*0.5*h;  
    F2f(EK2f,:)=F2f(EK2f,:)-[qsh*wt5,-qsh*wt5, qsh.*xib2*wt5]*0.5*h; 
    rZv=(rR+rRn)/2.*((zR-vR.^2/4-vRn.^2/4)-(zL-vL.^2/4-vLn.^2/4));
    F2f(EK2f,:)=F2f(EK2f,:)+[rZv*wt5,-rZv*wt5, rZv.*xib2*wt5]*0.5*h;
    qsh=-mu*(UL(8,:)+UL(9,:));
    F3f(EK1f,:)=F3f(EK1f,:)+[qsh*wt5, qsh*wt5, qsh.*xib1*wt5]*0.5*h;   
    F3f(EK2f,:)=F3f(EK2f,:)-[qsh*wt5,-qsh*wt5, qsh.*xib2*wt5]*0.5*h; 
    rZu=(rR+rRn)/8.*(uR+uRn).*(vR+vRn-vL-vLn);
    F3f(EK2f,:)=F3f(EK2f,:)+[rZu*wt5,-rZu*wt5, rZu.*xib2*wt5]*0.5*h; 
    qsh=UL(5,:);
    F4f(EK1f,:)=F4f(EK1f,:)+[qsh*wt5, qsh*wt5*3, qsh.*xib1*wt5*3]*0.5*h;   
    F4f(EK2f,:)=F4f(EK2f,:)-[qsh*wt5,-qsh*wt5*3, qsh.*xib2*wt5*3]*0.5*h; 
    qsh=(UR(1,:)+URn(1,:))/2;
    F5f(EK1f,:)=F5f(EK1f,:)+[qsh*wt5, qsh*wt5*3, qsh.*xib1*wt5*3]*0.5*h;   
    F5f(EK2f,:)=F5f(EK2f,:)-[qsh*wt5,-qsh*wt5*3, qsh.*xib2*wt5*3]*0.5*h; 
    qsh=(UR(2,:)+URn(2,:))/2;
    F7f(EK1f,:)=F7f(EK1f,:)+[qsh*wt5, qsh*wt5*3, qsh.*xib1*wt5*3]*0.5*h;   
    F7f(EK2f,:)=F7f(EK2f,:)-[qsh*wt5,-qsh*wt5*3, qsh.*xib2*wt5*3]*0.5*h; 
    qsh=(UR(3,:)+URn(3,:))/2;
    F9f(EK1f,:)=F9f(EK1f,:)+[qsh*wt5, qsh*wt5*3, qsh.*xib1*wt5*3]*0.5*h;   
    F9f(EK2f,:)=F9f(EK2f,:)-[qsh*wt5,-qsh*wt5*3, qsh.*xib2*wt5*3]*0.5*h; 
end

for k=1:fExn 
    id=fEx(k);   EK=Ex(id);   K1=Cells(EK.K1);    K2=Cells(EK.K2);   
    K1it=Cellsit(EK.K1);    K2it=Cellsit(EK.K2);   
    EK1f=K1.Indx_fine;  EK2f=K2.Indx_fine;  h=min(K1.Wid, K2.Wid);  
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
    f1h=1/4*(rT.*vT+rT.*vTn+rTn.*vT+rTn.*vTn);
    F1f(EK1f,:)=F1f(EK1f,:)+[f1h*wt5, f1h.*xib1*wt5*3, f1h*wt5*3]*0.5*h;  
    F1f(EK2f,:)=F1f(EK2f,:)-[f1h*wt5, f1h.*xib2*wt5*3,-f1h*wt5*3]*0.5*h; 
    qsh=-mu*(UB(8,:)+UB(9,:));
    F2f(EK1f,:)=F2f(EK1f,:)+[qsh*wt5, qsh.*xib1*wt5, qsh*wt5]*0.5*h;  
    F2f(EK2f,:)=F2f(EK2f,:)-[qsh*wt5, qsh.*xib2*wt5,-qsh*wt5]*0.5*h; 
    rZu=(rT+rTn)/8.*(vT+vTn).*(uT+uTn-uB-uBn);
    F2f(EK2f,:)=F2f(EK2f,:)+[rZu*wt5, rZu.*xib2*wt5,-rZu*wt5]*0.5*h; 
    qsh=-mu*(4/3*UB(10,:)-2/3*UB(7,:));
    F3f(EK1f,:)=F3f(EK1f,:)+[qsh*wt5, qsh.*xib1*wt5, qsh*wt5]*0.5*h;  
    F3f(EK2f,:)=F3f(EK2f,:)-[qsh*wt5, qsh.*xib2*wt5,-qsh*wt5]*0.5*h; 
    rZu=(rT+rTn)/2.*((zT-uT.^2/4-uTn.^2/4)-(zB-uB.^2/4-uBn.^2/4));
    F3f(EK2f,:)=F3f(EK2f,:)+[rZu*wt5, rZu.*xib2*wt5,-rZu*wt5]*0.5*h; 
    qsh=UB(6,:);
    F4f(EK1f,:)=F4f(EK1f,:)+[qsh*wt5, qsh.*xib1*wt5*3, qsh*wt5*3]*0.5*h;   
    F4f(EK2f,:)=F4f(EK2f,:)-[qsh*wt5, qsh.*xib2*wt5*3,-qsh*wt5*3]*0.5*h; 
    qsh=(UT(1,:)+UTn(1,:))/2;
    F6f(EK1f,:)=F6f(EK1f,:)+[qsh*wt5, qsh.*xib1*wt5*3, qsh*wt5*3]*0.5*h;   
    F6f(EK2f,:)=F6f(EK2f,:)-[qsh*wt5, qsh.*xib2*wt5*3,-qsh*wt5*3]*0.5*h; 
    qsh=(UT(2,:)+UTn(2,:))/2;
    F8f(EK1f,:)=F8f(EK1f,:)+[qsh*wt5, qsh.*xib1*wt5*3, qsh*wt5*3]*0.5*h;   
    F8f(EK2f,:)=F8f(EK2f,:)-[qsh*wt5, qsh.*xib2*wt5*3,-qsh*wt5*3]*0.5*h; 
    qsh=(UT(3,:)+UTn(3,:))/2;
    F10f(EK1f,:)=F10f(EK1f,:)+[qsh*wt5, qsh.*xib1*wt5*3, qsh*wt5*3]*0.5*h;   
    F10f(EK2f,:)=F10f(EK2f,:)-[qsh*wt5, qsh.*xib2*wt5*3,-qsh*wt5*3]*0.5*h; 
end
cuM=zeros(ko,ko);   Fs=zeros(fCellsn,nk);  Us=zeros(fCellsn, nk);
for k=1:fCellsn
    id=fCells(k);   Ckn=Cells(id);  h=Ckn.Wid;  Ukn=Ckn.U;  Uknxi=Ukn*gpM;
    rnxi=Uknxi(1,:);  unxi=Uknxi(2,:);   vnxi=Uknxi(3,:);   
    Uk=Cellsit(id).U;    Ukxi=Uk*gpM;  Zk=Uk(4,:);
    rxi=Ukxi(1,:);  uxi=Ukxi(2,:);   vxi=Ukxi(3,:);  ukt=Uk'; Us(k,:)=ukt(:);
      
    rkc=(rxi+rnxi)/2;    vkc=(vxi+vnxi)/2;    ukc=(uxi+unxi)/2;   
    for j=1:ko
        for i=1:ko
            cuM(i,j)=rkc.*gpM(i,:).*gpM(j,:)*wt/4;
        end
    end
    f1x=(rxi.*uxi+rnxi.*uxi+rxi.*unxi+rnxi.*unxi)/4;
    f1y=(rxi.*vxi+rnxi.*vxi+rxi.*vnxi+rnxi.*vnxi)/4;
    F1f(k,:)=F1f(k,:)-[0, f1x*wt, f1y*wt]*1.5*h;  F1f(k,:)=F1f(k,:)/h/h;
    f2x=-mu*(4/3*Ukxi(7,:)-2/3*Ukxi(10,:));  f2y=-mu*(Ukxi(8,:)+Ukxi(9,:));
    F2f(k,:)=F2f(k,:)-[0, f2x*wt, f2y*wt]*0.5*h;  %%
    f2x=rkc.*(Zk*tfx-vxi.*(Uk(3,:)*tfx)/2-vnxi.*(Ukn(3,:)*tfx)/2);  
    f2y=rkc.*vkc.*(Uk(2,:)*tfy+Ukn(2,:)*tfy)/2;
    F2f(k,:)=F2f(k,:)+[(f2x+f2y)*wt, (f2x+f2y).*gpM(2,:)*wt, (f2x+f2y).*gpM(3,:)*wt]*0.5*h;  
    F2f(k,:)=F2f(k,:)/h/h;
        
    f3y=-mu*(4/3*Ukxi(10,:)-2/3*Ukxi(7,:));  f3x=-mu*(Ukxi(8,:)+Ukxi(9,:));
    F3f(k,:)=F3f(k,:)-[0, f3x*wt, f3y*wt]*0.5*h;  
    f3y=rkc.*(Zk*tfy-uxi.*(Uk(2,:)*tfy)/2-unxi.*(Ukn(2,:)*tfy)/2);  
    f3x=rkc.*ukc.*(Uk(3,:)*tfx+Ukn(3,:)*tfx)/2;
    F3f(k,:)=F3f(k,:)+[(f3x+f3y)*wt, (f3x+f3y).*gpM(2,:)*wt, (f3x+f3y).*gpM(3,:)*wt]*0.5*h;  
    F3f(k,:)=F3f(k,:)/h/h;
      
    Fs(k,1:ko)=Uk(1,:)-Ukn(1,:)+dt*F1f(k,:);%%
    Fs(k,1+ko:2*ko)=(cuM*(Uk(2,:)'-Ukn(2,:)'))'+dt*F2f(k,:);%%
    Fs(k,1+2*ko:3*ko)=(cuM*(Uk(3,:)'-Ukn(3,:)'))'+dt*F3f(k,:);%%
    
    F4f(k,:)=F4f(k,:)-[0, Ukxi(5,:)*wt, Ukxi(6,:)*wt]*1.5*h;    F4f(k,:)=F4f(k,:)/h/h;
    zf=(unxi.^2+vnxi.^2+uxi.^2+vxi.^2)/4;%+Wr((rxi+rnxi)/2);
    if norm(Uk(1,:)-Ukn(1,:), inf)< myeps
        zf=zf+Wr((rxi+rnxi)/2);
    else
        zf=zf+(W(rxi)-W(rnxi))./(rxi-rnxi);
    end
    Fs(k,1+3*ko:4*ko)=Uk(4,:)-[zf*wt, zf.*gsx*wt*3, zf.*gsy*wt*3]/4+lambda*F4f(k,:);
    
    F5f(k,:)=F5f(k,:)-[0, (rxi+rnxi)/2*wt, 0]*1.5*h;   Fs(k,1+4*ko:5*ko)=Uk(5,:)-F5f(k,:)/h/h;
    F6f(k,:)=F6f(k,:)-[0, 0, (rxi+rnxi)/2*wt]*1.5*h;   Fs(k,1+5*ko:6*ko)=Uk(6,:)-F6f(k,:)/h/h;
    F7f(k,:)=F7f(k,:)-[0, (uxi+unxi)/2*wt, 0]*1.5*h;   Fs(k,1+6*ko:7*ko)=Uk(7,:)-F7f(k,:)/h/h;
    F8f(k,:)=F8f(k,:)-[0, 0, (uxi+unxi)/2*wt]*1.5*h;   Fs(k,1+7*ko:8*ko)=Uk(8,:)-F8f(k,:)/h/h;
    F9f(k,:)=F9f(k,:)-[0, (vxi+vnxi)/2*wt, 0]*1.5*h;   Fs(k,1+8*ko:9*ko)=Uk(9,:)-F9f(k,:)/h/h;
    F10f(k,:)=F10f(k,:)-[0, 0, (vxi+vnxi)/2*wt]*1.5*h; Fs(k,1+9*ko:10*ko)=Uk(10,:)-F10f(k,:)/h/h;
end
Fs=Fs';  Fu=Fs(:);
end