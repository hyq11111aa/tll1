function [Mesh, Mesh_fine]=Init(N,LEV,t)
global cell edge  mesh gs5 wt5 gsx gsy wt tf ko xRight lambda
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Initmesh %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MF=50;  ck=cell; 
Cells(1:N*N*MF)=cell;  Ex(1:N*N*MF)=edge;  Ey=Ex;  Eyn=0;  Exn=0;
%%%%%%%%%%%%%%%%%%%% initial mesh %%%%%%%%%%%%%%%%%%%%%%%%%%%
h=xRight/N;   Uk=zeros(10,ko);    Wes=sqrt(1/lambda);
% rvap=0.1;   rliq=0.6;
r1=0.25;  r2=0.1;  
rvap=0.3198;   rliq=1.8071;
% r1=0.15;   r2=0.05;  r3=0.07;   r4=0.03;  n=4;

A=load('xyrs'); xss=A(:,1); yss=A(:,2); rss=A(:,3);  n=50;
% n=50;
% xss=randperm(950000,n)*0.000001;   yss=randperm(950000,n)*0.000001;   
% rss=0.01+0.01*rand(1,n);
% fid=fopen('xyrs','wt');
% for k=1:n
%     fprintf(fid, '%8.6f %8.6f %8.6f \n', xss(k), yss(k), rss(k));
% end
% fclose(fid);
k=0;
for j=1:N
    for i=1:N
        k=k+1;
        %%%%%%%%%%%%%%%%%%%%% info of cells %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ck.Indx=k;  xo=(i-0.5)*h;   yo=(j-0.5)*h;   ck.Center=[xo, yo];
        ck.Wid=h;   xs=xo+0.5*h*gsx;    ys=yo+0.5*h*gsy;
%         d1=(xs-0.4).*(xs-0.4)+(ys-0.5).*(ys-0.5);     d1=sqrt(d1);
%         d2=(xs-0.78).*(xs-0.78)+(ys-0.5).*(ys-0.5);   d2=sqrt(d2);
%         r=0.1+0.25*(tanh(0.5*(d1-r1)*Wes)+tanh(0.5*(d2-r2)*Wes)); 
%         d1=sqrt((xs-0.3).^2+(ys-0.5).^2);     d2=sqrt((xs-0.7).^2+(ys-0.3).^2); 
%         d3=sqrt((xs-0.6).^2+(ys-0.7).^2);     d4=sqrt((xs-0.5).^2+(ys-0.5).^2); 
%         r=0.5*(n*rliq-(n-2)*rvap)+0.5*(rvap-rliq)*(tanh(0.5*(d1-r1)*Wes)+tanh(0.5*(d2-r2)*Wes)...
%             +tanh(0.5*(d3-r3)*Wes)+tanh(0.5*(d4-r4)*Wes));
        tans=0;
        for l=1:n
            dl=sqrt((xs-xss(l)).^2+(ys-yss(l)).^2); 
            tans=tans+tanh(0.5*(dl-rss(l))*Wes);
        end    
        r=0.5*(n*rliq-(n-2)*rvap)+0.5*(rvap-rliq)*tans;  
        Uk(1,:)=[r*wt, r.*gsx*wt*3, r.*gsy*wt*3]/4; ck.U=Uk;
        Cells(k)=ck;
        %%%%%%%%%%%%%%%%%%%% edges in x direction %%%%%%%%%%%%%%%%%%%%%%%%
        Exn=Exn+1;  ek=edge;   ek.Indx=Exn;  ek.K1=(j-2)*N+i;  ek.K2=(j-1)*N+i;     
        if j==1
            ek.K1=(N-1)*N+i; 
        end
        Ex(Exn)=ek;
        
        %%%%%%%%%%%%%%%%%%%% edges in y direction %%%%%%%%%%%%%%%%%%%%%%%%
        Eyn=Eyn+1;  ek=edge;   ek.Indx=Eyn;   ek.K1=(j-1)*N+i-1; ek.K2=(j-1)*N+i; 
        if i==1
            ek.K1=j*N;
        end
        Ey(Eyn)=ek;
    end
end
for k=1:Exn
    ex=Ex(k);       Cells(ex.K1).Edges(4)=ex.Indx;   Cells(ex.K2).Edges(3)=ex.Indx;
end
for k=1:Eyn
    ey=Ey(k);       Cells(ey.K1).Edges(2)=ey.Indx;   Cells(ey.K2).Edges(1)=ey.Indx;
end
Mesh=mesh;
Mesh.Ex=Ex;  Mesh.Exn=Exn;   Mesh.Ey=Ey;  Mesh.Eyn=Eyn;
Mesh.Cells=Cells; Mesh.Cellsn=N*N;
[Mesh, Mesh_fine]=SettingC(Mesh);
% %%%%%%%%%%%%%%%%%%%%%%% finest cells, edges %%%%%%%%%%%%%
for l=1:LEV
    [Mesh, TroubCells, Troubn]=Trouble_Cellx(Mesh, Mesh_fine, LEV); 
    disp(['Troubn=' num2str(Troubn)])
    [Mesh]=RefineT0(Mesh, TroubCells, Troubn);
    [Mesh, Mesh_fine]=SettingC(Mesh);
end