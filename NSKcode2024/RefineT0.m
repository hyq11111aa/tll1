function   [Mesh]=RefineT0(Mesh, TroubCells, Troubn)
global cell edge gsx gsy wt lambda
Cells=Mesh.Cells;    Cellsn=Mesh.Cellsn;    
Ex=Mesh.Ex;   Exn=Mesh.Exn; Ey=Mesh.Ey;    Eyn=Mesh.Eyn;  Wes=sqrt(1/lambda);  
% r1=0.25;   r2=0.1;
% rvap=0.1;   rliq=0.6;
rvap=0.3198;   rliq=1.8071;
% r1=0.15;   r2=0.05;  r3=0.07;   r4=0.03;  n=4; 
A=load('xyrs'); xss=A(:,1); yss=A(:,2); rss=A(:,3); n=50;  
%%%%%%%%%%%%%%%%%%%%%%%%% refine trouble cells %%%%%%%%%%%%%%%%%%%%%%%%%%%
for k=1:Troubn
    Cell=Cells(TroubCells(k));   h=Cell.Wid;   center=Cell.Center;
    Cell.HaveChild=1;       Clev=Cell.Lev;  
    
    %%%%%%%%%%%%%%%%%%%%%%%% divide a trouble element%%%%%%%%%%%%%%%%%%%%%
    c1=cell;    c1.Indx=Cellsn+1;    c1.Lev=Clev+1;    c1.Wid=h/2;   c1.Parent=Cell.Indx;   
    xo=center(1)-h/4;    yo=center(2)-h/4;  c1.Center=[xo, yo];
    xv=xo+0.5*h/2*gsx;    yv=yo+0.5*h/2*gsy;
%     d1=(xv-0.4).*(xv-0.4)+(yv-0.5).*(yv-0.5);     d1=sqrt(d1);
%     d2=(xv-0.78).*(xv-0.78)+(yv-0.5).*(yv-0.5);   d2=sqrt(d2);
%     r=0.1+0.25*(tanh(0.5*(d1-r1)*Wes)+tanh(0.5*(d2-r2)*Wes)); 
%     d1=sqrt((xv-0.3).^2+(yv-0.5).^2);     d2=sqrt((xv-0.7).^2+(yv-0.3).^2); 
%     d3=sqrt((xv-0.6).^2+(yv-0.7).^2);     d4=sqrt((xv-0.5).^2+(yv-0.5).^2); 
%     r=0.5*(n*rliq-(n-2)*rvap)+0.5*(rvap-rliq)*(tanh(0.5*(d1-r1)*Wes)+tanh(0.5*(d2-r2)*Wes)...
%             +tanh(0.5*(d3-r3)*Wes)+tanh(0.5*(d4-r4)*Wes));
    tans=0;
    for l=1:n
        dl=sqrt((xv-xss(l)).^2+(yv-yss(l)).^2); 
        tans=tans+tanh(0.5*(dl-rss(l))*Wes);
    end    
    r=0.5*(n*rliq-(n-2)*rvap)+0.5*(rvap-rliq)*tans;   
    c1.U(1,:)=[r*wt, r.*gsx*wt*3, r.*gsy*wt*3]/4;
   
    c2=c1;     c2.Indx=Cellsn+2;  
    xo=center(1)+h/4;    yo=center(2)-h/4;  c2.Center=[xo, yo];
    xs=xo+0.5*h/2*gsx;    ys=yo+0.5*h/2*gsy;
%     d1=(xs-0.4).*(xs-0.4)+(ys-0.5).*(ys-0.5);     d1=sqrt(d1);
%     d2=(xs-0.78).*(xs-0.78)+(ys-0.5).*(ys-0.5);   d2=sqrt(d2);
%     r=0.1+0.25*(tanh(0.5*(d1-r1)*Wes)+tanh(0.5*(d2-r2)*Wes)); 
%     d1=sqrt((xs-0.3).^2+(ys-0.5).^2);     d2=sqrt((xs-0.7).^2+(ys-0.3).^2); 
%     d3=sqrt((xs-0.6).^2+(ys-0.7).^2);     d4=sqrt((xs-0.5).^2+(ys-0.5).^2); 
%     r=0.5*(n*rliq-(n-2)*rvap)+0.5*(rvap-rliq)*(tanh(0.5*(d1-r1)*Wes)+tanh(0.5*(d2-r2)*Wes)...
%             +tanh(0.5*(d3-r3)*Wes)+tanh(0.5*(d4-r4)*Wes));
    tans=0;
    for l=1:n
        dl=sqrt((xs-xss(l)).^2+(ys-yss(l)).^2); 
        tans=tans+tanh(0.5*(dl-rss(l))*Wes);
    end    
    r=0.5*(n*rliq-(n-2)*rvap)+0.5*(rvap-rliq)*tans;  
    c2.U(1,:)=[r*wt, r.*gsx*wt*3, r.*gsy*wt*3]/4;
    
    c3=c1;     c3.Indx=Cellsn+3; 
    xo=center(1)+h/4;    yo=center(2)+h/4;  c3.Center=[xo, yo];
    xs=xo+0.5*h/2*gsx;    ys=yo+0.5*h/2*gsy;
%     d1=(xs-0.4).*(xs-0.4)+(ys-0.5).*(ys-0.5);     d1=sqrt(d1);
%     d2=(xs-0.78).*(xs-0.78)+(ys-0.5).*(ys-0.5);   d2=sqrt(d2);
%     r=0.1+0.25*(tanh(0.5*(d1-r1)*Wes)+tanh(0.5*(d2-r2)*Wes)); 
%     d1=sqrt((xs-0.3).^2+(ys-0.5).^2);     d2=sqrt((xs-0.7).^2+(ys-0.3).^2); 
%     d3=sqrt((xs-0.6).^2+(ys-0.7).^2);     d4=sqrt((xs-0.5).^2+(ys-0.5).^2); 
%     r=0.5*(n*rliq-(n-2)*rvap)+0.5*(rvap-rliq)*(tanh(0.5*(d1-r1)*Wes)+tanh(0.5*(d2-r2)*Wes)...
%             +tanh(0.5*(d3-r3)*Wes)+tanh(0.5*(d4-r4)*Wes));
    tans=0;
    for l=1:n
        dl=sqrt((xs-xss(l)).^2+(ys-yss(l)).^2); 
        tans=tans+tanh(0.5*(dl-rss(l))*Wes);
    end    
    r=0.5*(n*rliq-(n-2)*rvap)+0.5*(rvap-rliq)*tans;  
    c3.U(1,:)=[r*wt, r.*gsx*wt*3, r.*gsy*wt*3]/4;
    
    c4=c1;     c4.Indx=Cellsn+4;  
    xo=center(1)-h/4;    yo=center(2)+h/4;  c4.Center=[xo, yo];
    xs=xo+0.5*h/2*gsx;    ys=yo+0.5*h/2*gsy;
%     d1=(xs-0.4).*(xs-0.4)+(ys-0.5).*(ys-0.5);     d1=sqrt(d1);
%     d2=(xs-0.78).*(xs-0.78)+(ys-0.5).*(ys-0.5);   d2=sqrt(d2);
%     r=0.1+0.25*(tanh(0.5*(d1-r1)*Wes)+tanh(0.5*(d2-r2)*Wes)); 
%     d1=sqrt((xs-0.3).^2+(ys-0.5).^2);     d2=sqrt((xs-0.7).^2+(ys-0.3).^2); 
%     d3=sqrt((xs-0.6).^2+(ys-0.7).^2);     d4=sqrt((xs-0.5).^2+(ys-0.5).^2); 
%     r=0.5*(n*rliq-(n-2)*rvap)+0.5*(rvap-rliq)*(tanh(0.5*(d1-r1)*Wes)+tanh(0.5*(d2-r2)*Wes)...
%             +tanh(0.5*(d3-r3)*Wes)+tanh(0.5*(d4-r4)*Wes));
    tans=0;
    for l=1:n
        dl=sqrt((xs-xss(l)).^2+(ys-yss(l)).^2); 
        tans=tans+tanh(0.5*(dl-rss(l))*Wes);
    end    
    r=0.5*(n*rliq-(n-2)*rvap)+0.5*(rvap-rliq)*tans;  
    c4.U(1,:)=[r*wt, r.*gsx*wt*3, r.*gsy*wt*3]/4;

    %%%%%%%%%%%%%%%%%%%% decide whether a node exist %%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%% if not exist, add node and new edges %%%%%%%%%%%%%%%
    eds=Cell.Edges; 
    eyL=Ey(eds(1));  eyR=Ey(eds(2));  exB=Ex(eds(3));  exT=Ex(eds(4));
    if ~exB.HaveChild
        iE=eds(3);    Ex1=edge;   Ex2=edge;   
        Ex1.Indx=Exn+1;      Ex1.Parent=iE;   Ex1.Lev=Clev+1;    
        c1.Edges(3)=Exn+1;   Ex1.K1=exB.K1;   Ex1.K2=c1.Indx;
        Ex2.Indx=Exn+2;      Ex2.Parent=iE;   Ex2.Lev=Clev+1;    
        c2.Edges(3)=Exn+2;   Ex2.K1=exB.K1;   Ex2.K2=c2.Indx;
        Ex(iE).HaveChild=1;  Ex(iE).Cs=[Exn+1,Exn+2];   
        Ex(Exn+1)=Ex1;       Ex(Exn+2)=Ex2;   Exn=Exn+2;
    else
        exBc=exB.Cs;    exBc1=exBc(1);    exBc2=exBc(2);
        c1.Edges(3)=exBc1;  c2.Edges(3)=exBc2;
        Ex(exBc1).K2=c1.Indx;   Ex(exBc2).K2=c2.Indx;
    end
    %%%%%%%% left boundary
    if ~eyL.HaveChild
        iE=Cell.Edges(1);   Ey1=edge;   Ey2=edge;
        Ey1.Indx=Eyn+1;     Ey1.Parent=iE;  Ey1.Lev=Clev+1;   
        c1.Edges(1)=Eyn+1;  Ey1.K1=eyL.K1;  Ey1.K2=c1.Indx;
        Ey2.Indx=Eyn+2;     Ey2.Parent=iE;  Ey2.Lev=Clev+1;   
        c4.Edges(1)=Eyn+2;  Ey2.K1=eyL.K1;  Ey2.K2=c4.Indx;
        Ey(iE).HaveChild=1;     Ey(iE).Cs=[Eyn+1,Eyn+2];
        Ey(Eyn+1)=Ey1;   Ey(Eyn+2)=Ey2;   Eyn=Eyn+2;
    else
        eyLc=eyL.Cs;    eyLc1=eyLc(1);    eyLc2=eyLc(2);
        c1.Edges(1)=eyLc1;  c4.Edges(1)=eyLc2;
        Ey(eyLc1).K2=c1.Indx;   Ey(eyLc2).K2=c4.Indx;
    end
    
    Ey1=edge;   Ey2=edge; Ex1=edge;   Ex2=edge;  
    Ex1.Indx=Exn+1;      Ex1.K1=c1.Indx;   Ex1.K2=c4.Indx;   
    Ex1.Lev=Clev+1;      c1.Edges(4)=Exn+1;   c4.Edges(3)=Exn+1;
    Ex2.Indx=Exn+2;      Ex2.K1=c2.Indx;   Ex2.K2=c3.Indx;     
    Ex2.Lev=Clev+1;      c2.Edges(4)=Exn+2;   c3.Edges(3)=Exn+2;
    Ex(Exn+1)=Ex1;       Ex(Exn+2)=Ex2;      Exn=Exn+2;     
           
    Ey1.Indx=Eyn+1;      Ey1.K1=c1.Indx;     Ey1.K2=c2.Indx;     
    Ey1.Lev=Clev+1;      c1.Edges(2)=Eyn+1;  c2.Edges(1)=Eyn+1;   
    Ey2.Indx=Eyn+2;      Ey2.K1=c4.Indx;     Ey2.K2=c3.Indx;       
    Ey2.Lev=Clev+1;      c4.Edges(2)=Eyn+2;  c3.Edges(1)=Eyn+2;  
    Ey(Eyn+1)=Ey1;       Ey(Eyn+2)=Ey2;      Eyn=Eyn+2;  
    
     if ~eyR.HaveChild
        iE=eds(2);    Ey1=edge;   Ey2=edge;  
        Ey1.Indx=Eyn+1;   Ey1.Parent=iE;    Ey1.Lev=Clev+1;
        c2.Edges(2)=Eyn+1;  Ey1.K2=eyR.K2;  Ey1.K1=c2.Indx;
        Ey2.Indx=Eyn+2;   Ey2.Parent=iE;    Ey2.Lev=Clev+1;
        c3.Edges(2)=Eyn+2;  Ey2.K2=eyR.K2;  Ey2.K1=c3.Indx;
        Ey(iE).HaveChild=1;  Ey(iE).Cs=[Eyn+1,Eyn+2];
        Ey(Eyn+1)=Ey1;       Ey(Eyn+2)=Ey2;   Eyn=Eyn+2;
     else
        eyRc=eyR.Cs;    eyc1=eyRc(1);    eyc2=eyRc(2);
        c2.Edges(2)=eyc1;  c3.Edges(2)=eyc2;
        Ey(eyc1).K1=c2.Indx;   Ey(eyc2).K1=c3.Indx;
     end
     %%%%%%%% right boundary
     if ~exT.HaveChild
        iE=eds(4);    Ex1=edge;   Ex2=edge; 
        Ex1.Indx=Exn+1;      Ex1.Parent=iE;     Ex1.Lev=Clev+1;
        c4.Edges(4)=Exn+1;   Ex1.K2=exT.K2;     Ex1.K1=c4.Indx;       
        Ex2.Indx=Exn+2;      Ex2.Parent=iE;     Ex2.Lev=Clev+1;
        c3.Edges(4)=Exn+2;   Ex2.K2=exT.K2;     Ex2.K1=c3.Indx;
        Ex(iE).HaveChild=1;  Ex(iE).Cs=[Exn+1,Exn+2];
        Ex(Exn+1)=Ex1;       Ex(Exn+2)=Ex2;   Exn=Exn+2;
     else
        exTc=exT.Cs;    exc1=exTc(1);    exc2=exTc(2);
        c4.Edges(4)=exc1;  c3.Edges(4)=exc2;
        Ex(exc1).K1=c4.Indx;   Ex(exc2).K1=c3.Indx;
     end
    
    Cell.Cs=[Cellsn+1, Cellsn+2, Cellsn+3, Cellsn+4];  Cells(Cell.Indx)=Cell;
    Cells(Cellsn+1)=c1;       Cells(Cellsn+2)=c2;
    Cells(Cellsn+3)=c3;       Cells(Cellsn+4)=c4;
    Cellsn=Cellsn+4;
end

for k=1:Cellsn
     Cells(k).NoTrouble=1;
end
Mesh.Ex=Ex;  Mesh.Exn=Exn;   Mesh.Ey=Ey;  Mesh.Eyn=Eyn;
Mesh.Cells=Cells;   Mesh.Cellsn=Cellsn;
    