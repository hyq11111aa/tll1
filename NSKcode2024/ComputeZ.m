function [Cellsit]=ComputeZ(Cellsit, Mesh, Mesh_fine)
global gs5 wt5 gps1 gps2 wt gsx gsy ko lambda Wr
global gpL gpR gpB gpT  gcL1 gcR1 gcB1 gcT1  gcL2 gcR2 gcB2 gcT2 gpM
%%%%%%%%%%%%%%%%%%%%% Mesh_fine u1,u2,p
fCells=Mesh_fine.fCells;    fCellsn=Mesh_fine.fCellsn;
fEy=Mesh_fine.fEy;   fEyn=Mesh_fine.fEyn;   
fEx=Mesh_fine.fEx;   fExn=Mesh_fine.fExn;             
Cells=Mesh.Cells;   Ey=Mesh.Ey; Ex=Mesh.Ex;

q1f=zeros(fCellsn,ko);  q2f=q1f;  q3f=q1f;  q4f=q1f;  q5f=q1f;  q6f=q1f;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% q_1=rho_x, q_2=rho_y, q_3=u_x, q_4=u_y, q_5=v_x, q_6=v_y
for k=1:fEyn 
    id=fEy(k);   EK=Ey(id);   K1=Cells(EK.K1);    K2=Cells(EK.K2);   K2it=Cellsit(EK.K2);
    EK1f=K1.Indx_fine;  EK2f=K2.Indx_fine;  h=min(K1.Wid, K2.Wid); 
    if K1.Lev==K2.Lev
         xib1=gs5;    UR=K2.U*gpR;  URit=K2it.U*gpR;  xib2=gs5;
    elseif K1.Lev==K2.Lev-1
        if K1.Center(2)>K2.Center(2)
            xib1=gps1; 
        else
            xib1=gps2;
        end
        UR=K2.U*gpR;  URit=K2it.U*gpR;  xib2=gs5;
    elseif K1.Lev==K2.Lev+1
        if K2.Center(2)>K1.Center(2)
            UR=K2.U*gcR1;  xib2=gps1;   URit=K2it.U*gcR1;
        else
            UR=K2.U*gcR2;  xib2=gps2;   URit=K2it.U*gcR2;
        end
        xib1=gs5;  
    else
        disp('Init: wrong conforming')
    end  
    rh=(UR(1,:)+URit(1,:))/2;  uh=(UR(2,:)+URit(2,:))/2;  vh=(UR(3,:)+URit(3,:))/2;  
    q1f(EK1f,:)=q1f(EK1f,:)+[rh*wt5, rh*wt5*3, rh.*xib1*wt5*3]*0.5*h;  
    q1f(EK2f,:)=q1f(EK2f,:)-[rh*wt5,-rh*wt5*3, rh.*xib2*wt5*3]*0.5*h; 
    q3f(EK1f,:)=q3f(EK1f,:)+[uh*wt5, uh*wt5*3, uh.*xib1*wt5*3]*0.5*h;  
    q3f(EK2f,:)=q3f(EK2f,:)-[uh*wt5,-uh*wt5*3, uh.*xib2*wt5*3]*0.5*h; 
    q5f(EK1f,:)=q5f(EK1f,:)+[vh*wt5, vh*wt5*3, vh.*xib1*wt5*3]*0.5*h;  
    q5f(EK2f,:)=q5f(EK2f,:)-[vh*wt5,-vh*wt5*3, vh.*xib2*wt5*3]*0.5*h; 
end

for k=1:fExn 
    id=fEx(k);   EK=Ex(id);   K1=Cells(EK.K1);    K2=Cells(EK.K2);   K2it=Cellsit(EK.K2);   
    EK1f=K1.Indx_fine;  EK2f=K2.Indx_fine;  h=min(K1.Wid, K2.Wid);  
    if K1.Lev==K2.Lev
        xib1=gs5;   UT=K2.U*gpT;   UTit=K2it.U*gpT;  xib2=gs5; 
    elseif K1.Lev==K2.Lev-1
        if K1.Center(1)>K2.Center(1)
            xib1=gps1;    
        else
            xib1=gps2;   
        end
        UT=K2.U*gpT;   UTit=K2it.U*gpT;    xib2=gs5;
    elseif K1.Lev==K2.Lev+1
        if K2.Center(1)>K1.Center(1)
           UT=K2.U*gcT1;   xib2=gps1;   UTit=K2it.U*gcT1; 
        else
           UT=K2.U*gcT2;   xib2=gps2;   UTit=K2it.U*gcT2; 
        end
        xib1=gs5;
    else
        disp('Init: wrong conforming')
    end 
    rh=(UT(1,:)+UTit(1,:))/2; uh=(UT(2,:)+UTit(2,:))/2; vh=(UT(3,:)+UTit(3,:))/2; 
    q2f(EK1f,:)=q2f(EK1f,:)+[rh*wt5, rh.*xib1*wt5*3, rh*wt5*3]*0.5*h;  
    q2f(EK2f,:)=q2f(EK2f,:)-[rh*wt5, rh.*xib2*wt5*3,-rh*wt5*3]*0.5*h; 
    q4f(EK1f,:)=q4f(EK1f,:)+[uh*wt5, uh.*xib1*wt5*3, uh*wt5*3]*0.5*h;  
    q4f(EK2f,:)=q4f(EK2f,:)-[uh*wt5, uh.*xib2*wt5*3,-uh*wt5*3]*0.5*h; 
    q6f(EK1f,:)=q6f(EK1f,:)+[vh*wt5, vh.*xib1*wt5*3, vh*wt5*3]*0.5*h;  
    q6f(EK2f,:)=q6f(EK2f,:)-[vh*wt5, vh.*xib2*wt5*3,-vh*wt5*3]*0.5*h; 
end
for k=1:fCellsn
    id=fCells(k);   Ck=Cells(id); Ckit=Cellsit(id); h=Ck.Wid;  
    Ukxi=Ck.U*gpM;  Ukitxi=Ckit.U*gpM;   rxi=(Ukxi(1,:)+Ukitxi(1,:))/2;  
    q1f(k,:)=q1f(k,:)-[0, rxi*wt, 0]*1.5*h;  q1f(k,:)=q1f(k,:)/h/h;
    q2f(k,:)=q2f(k,:)-[0, 0, rxi*wt]*1.5*h;  q2f(k,:)=q2f(k,:)/h/h;
    uxi=(Ukxi(2,:)+Ukitxi(2,:))/2;  
    q3f(k,:)=q3f(k,:)-[0, uxi*wt, 0]*1.5*h;  q3f(k,:)=q3f(k,:)/h/h;
    q4f(k,:)=q4f(k,:)-[0, 0, uxi*wt]*1.5*h;  q4f(k,:)=q4f(k,:)/h/h;
    vxi=(Ukxi(3,:)+Ukitxi(3,:))/2;  
    q5f(k,:)=q5f(k,:)-[0, vxi*wt, 0]*1.5*h;  q5f(k,:)=q5f(k,:)/h/h;
    q6f(k,:)=q6f(k,:)-[0, 0, vxi*wt]*1.5*h;  q6f(k,:)=q6f(k,:)/h/h;
    Cellsit(id).U(5,:)=q1f(k,:);    Cellsit(id).U(6,:)=q2f(k,:);
    Cellsit(id).U(7,:)=q3f(k,:);    Cellsit(id).U(8,:)=q4f(k,:);
    Cellsit(id).U(9,:)=q5f(k,:);    Cellsit(id).U(10,:)=q6f(k,:);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Z-lambda(r_xx+r_yy)=Z-lambda(q_1x+q_2y)
F1f=zeros(fCellsn,ko);
for k=1:fEyn 
    id=fEy(k);   EK=Ey(id);   K1=Cellsit(EK.K1);    K2=Cellsit(EK.K2);      
    EK1f=K1.Indx_fine;  EK2f=K2.Indx_fine;  h=min(K1.Wid, K2.Wid);  
    if K1.Lev==K2.Lev
         xib1=gs5;    q1L=K1.U(5,:)*gpL;   xib2=gs5;
    elseif K1.Lev==K2.Lev-1
        if K1.Center(2)>K2.Center(2)
            xib1=gps1;  q1L=K1.U(5,:)*gcL1; 
        else
            xib1=gps2;  q1L=K1.U(5,:)*gcL2; 
        end
        xib2=gs5;
    elseif K1.Lev==K2.Lev+1
        if K2.Center(2)>K1.Center(2)
            xib2=gps1;   
        else
            xib2=gps2;    
        end
        xib1=gs5;  q1L=K1.U(5,:)*gpL; 
    else
        disp('Init: wrong conforming')
    end  
    F1f(EK1f,:)=F1f(EK1f,:)+[q1L*wt5, q1L*wt5*3, q1L.*xib1*wt5*3]*0.5*h;  
    F1f(EK2f,:)=F1f(EK2f,:)-[q1L*wt5,-q1L*wt5*3, q1L.*xib2*wt5*3]*0.5*h; 
end

for k=1:fExn 
    id=fEx(k);   EK=Ex(id);   K1=Cellsit(EK.K1);    K2=Cellsit(EK.K2);   
    EK1f=K1.Indx_fine;  EK2f=K2.Indx_fine;  h=min(K1.Wid, K2.Wid);  
    if K1.Lev==K2.Lev
        xib1=gs5;   q2B=K1.U(6,:)*gpB;  xib2=gs5; 
    elseif K1.Lev==K2.Lev-1
        if K1.Center(1)>K2.Center(1)
            xib1=gps1;    q2B=K1.U(6,:)*gcB1;
        else
            xib1=gps2;    q2B=K1.U(6,:)*gcB2;
        end
        xib2=gs5;
    elseif K1.Lev==K2.Lev+1
        if K2.Center(1)>K1.Center(1)
           xib2=gps1;   
        else
           xib2=gps2;   
        end
        xib1=gs5;   q2B=K1.U(6,:)*gpB;
    else
        disp('Init: wrong conforming')
    end 
    F1f(EK1f,:)=F1f(EK1f,:)+[q2B*wt5, q2B.*xib1*wt5*3, q2B*wt5*3]*0.5*h;  
    F1f(EK2f,:)=F1f(EK2f,:)-[q2B*wt5, q2B.*xib2*wt5*3,-q2B*wt5*3]*0.5*h; 
end
for k=1:fCellsn
    id=fCells(k);   Ck=Cells(id);  Ckit=Cellsit(id);  h=Ck.Wid;   Ukit=Ckit.U;  
    q1xi=Ukit(5,:)*gpM;    q2xi=Ukit(6,:)*gpM; 
    F1f(k,:)=F1f(k,:)-[0, q1xi*wt, q2xi*wt]*1.5*h;  F1f(k,:)=F1f(k,:)/h/h;
    Ck_Uxi=Ck.U*gpM;  rxin=Ck_Uxi(1,:);  uxin=Ck_Uxi(2,:);   vxin=Ck_Uxi(3,:);
    rxi=Ukit(1,:)*gpM;    uxi=Ukit(2,:)*gpM;    vxi=Ukit(3,:)*gpM;
    zf=(uxi.^2+vxi.^2+uxin.^2+vxin.^2)/4+Wr((rxi+rxin)/2);
    Cellsit(id).U(4,:)=[zf*wt; zf.*gsx*wt*3; zf.*gsy*wt*3]/4-lambda*F1f(k,:)';
end
     
end