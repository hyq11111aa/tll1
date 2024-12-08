function  MEs=MEfun(Mesh, Mesh_fine, t)
global  gpM wt W lambda
Mass=0;  Energy=0;  
fCells=Mesh_fine.fCells;   fCellsn=Mesh_fine.fCellsn;   Cells=Mesh.Cells;    
for k=1:fCellsn
    idk=fCells(k);   Ck=Cells(idk);   Ukxi=Ck.U*gpM;    h=Ck.Wid;
    rxi=Ukxi(1,:);   uxi=Ukxi(2,:);   vxi=Ukxi(3,:);  q1xi=Ukxi(5,:);  q2xi=Ukxi(6,:);
    Mass=Mass+(rxi*wt*h*h)/4;  
    Energy=Energy+((W(rxi)+1/2*lambda*(q1xi.*q1xi+q2xi.*q2xi)+0.5*rxi.*(uxi.*uxi+vxi.*vxi))*wt*h*h)/4;
end
MEs=[t, Mass, Energy, fCellsn];
end 