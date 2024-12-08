function  [Mesh, Mesh_fine, MEv]=LDG_scheme(Mesh, Mesh_fine, LEV, T, N)
global xRight  gpM wt W lambda
h=xRight/N;
t=0;   dt=0.256*h/2^LEV;  it=1;
MEv=zeros(10^8,4);  
Mass=0;  Energy=0;  
fCells=Mesh_fine.fCells;   fCellsn=Mesh_fine.fCellsn;   Cells=Mesh.Cells;    
for k=1:fCellsn
    idk=fCells(k);   Ck=Cells(idk);   Ukxi=Ck.U*gpM;    h=Ck.Wid;
    rxi=Ukxi(1,:);   uxi=Ukxi(2,:);   vxi=Ukxi(3,:);  q1xi=Ukxi(5,:);  q2xi=Ukxi(6,:);
    Mass=Mass+(rxi*wt*h*h)/4;  
    Energy=Energy+((W(rxi)+1/2*lambda*(q1xi.*q1xi+q2xi.*q2xi)+0.5*rxi.*(uxi.*uxi+vxi.*vxi))*wt*h*h)/4;
end
MEv(it,:)=[t, Mass, Energy, fCellsn];

while t<T
    if t+dt>=T
        dt=T-t;
    end
    if rem(it,10)==0
        [Mesh, TroubCells, Troubn]=Trouble_Cellx(Mesh, Mesh_fine, LEV);
        Troubn
        if Troubn>0
            [Mesh]=Refine(Mesh, TroubCells, Troubn);
            [Mesh, Mesh_fine]=SettingC(Mesh);
        end
        [CoarsCs, Coarsnum]=CTrouble_Cell(Mesh, Mesh_fine);
        Coarsnum
        if Coarsnum>0
            [Mesh]=Coarsenx(Mesh, Coarsnum, CoarsCs);
            [Mesh, Mesh_fine]=SettingC(Mesh);
        end
    end
    
     if rem(it,500)==0
        Cells=Mesh.Cells;  fCells=Mesh_fine.fCells;  fCellsn=Mesh_fine.fCellsn;
        filenam=['B2fullyL' num2str(LEV) 't' num2str(t) '.dat'];
        Savedata(Cells, fCells, fCellsn, filenam);
     end
    [Mesh]=Time2Newton(Mesh, Mesh_fine, t, dt);
        
    t=t+dt;  it=it+1; disp(['t=' num2str(t)] )
    
    Mass=0;  Energy=0;  
    fCells=Mesh_fine.fCells;   fCellsn=Mesh_fine.fCellsn;   Cells=Mesh.Cells;    
    for k=1:fCellsn
        idk=fCells(k);   Ck=Cells(idk);   Ukxi=Ck.U*gpM;    h=Ck.Wid;
        rxi=Ukxi(1,:);   uxi=Ukxi(2,:);   vxi=Ukxi(3,:);  q1xi=Ukxi(5,:);  q2xi=Ukxi(6,:);
        Mass=Mass+(rxi*wt*h*h)/4;  
        Energy=Energy+((W(rxi)+1/2*lambda*(q1xi.*q1xi+q2xi.*q2xi)+0.5*rxi.*(uxi.*uxi+vxi.*vxi))*wt*h*h)/4;
    end
    MEv(it,:)=[t, Mass, Energy, fCellsn];
end
MEv=MEv(1:it,:);