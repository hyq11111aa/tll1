function [Mesh, TroubCells, Troubn]=Trouble_Cellx(Mesh, Mesh_fine, LEV)

fCells=Mesh_fine.fCells;    fCellsn=Mesh_fine.fCellsn;  
Cells=Mesh.Cells;   Ey=Mesh.Ey;        Ex=Mesh.Ex;     
%%%%%%%%%%%%%%%%%%%%%%% look for trouble cells %%%%%%%%%%%%%%%%%%%%%%%%%%%
eta=zeros(1,fCellsn);  
for k=1:fCellsn
    Ck=Cells(fCells(k));   
    rk=Ck.U(1,:);  etac=sqrt(2.0)*2.0*sqrt(rk(2)*rk(2)+rk(3)*rk(3));   
    eta(k)=etac;    
end
eta_new=eta;
%%%%%%%%%%%%%%%%%%%%%% m=2 %%%%%%%%%%%%%%%%%%%%%%%%
for k=1:fCellsn
  id=fCells(k);    Cc=Cells(id);    em=eta(k);      
  for i=1:Cc.Neighb_num                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     
       xc=Cc.Neighb(i);   Ccn=Cells(xc);  
	       em=max(eta(Ccn.Indx_fine),em);    
       for j=1:Ccn.Neighb_num
           xc=Ccn.Neighb(j);     Ccnx=Cells(xc);    
          em=max(eta(Ccnx.Indx_fine),em);	  
       end
  end
  eta_new(k)=em;   Cells(id).eta=em;
end
eta=eta_new;

eta1=0.05;%%1.0;%
Troubn=0;  TroubCells=zeros(1,floor(fCellsn/2));
for k=1:fCellsn
    fk=fCells(k);
    if eta(k)>eta1 && (Cells(fk).Lev<LEV)
        Troubn=Troubn+1;  Cells(fk).NoTrouble=0;
        TroubCells(Troubn)=fk; 
    end
end

%%%%%%%%%%%%%%%%%%%%%%Refinement: make sure conform level of mesh %%%%%%%%%
k=1;
while k<=Troubn
    Cell=Cells(TroubCells(k));     
    EL=Cell.Edges(1);  nKL=Ey(EL).K1;   CellL=Cells(nKL); 
    if (Cell.Lev-CellL.Lev>0) && CellL.NoTrouble
        Troubn=Troubn+1;  Cells(nKL).NoTrouble=0;   TroubCells(Troubn)=nKL;
    end
    ER=Cell.Edges(2);  nKR=Ey(ER).K2;   CellR=Cells(nKR); 
    if (Cell.Lev-CellR.Lev>0) && CellR.NoTrouble
        Troubn=Troubn+1;  Cells(nKR).NoTrouble=0;    TroubCells(Troubn)=nKR;
    end
    EB=Cell.Edges(3);  nKB=Ex(EB).K1;  CellB=Cells(nKB);  
    if (Cell.Lev-CellB.Lev>0) && CellB.NoTrouble
        Troubn=Troubn+1;  Cells(nKB).NoTrouble=0;   TroubCells(Troubn)=nKB;
    end
    ET=Cell.Edges(4);   nKT=Ex(ET).K2;  CellT=Cells(nKT);
    if (Cell.Lev-CellT.Lev>0) && CellT.NoTrouble
       Troubn=Troubn+1;  Cells(nKT).NoTrouble=0;  TroubCells(Troubn)=nKT;
    end
        
    k=k+1;
end
Mesh.Cells=Cells;

for i=1:Troubn
    for j=i+1:Troubn
        if TroubCells(j)<TroubCells(i)
            amin=TroubCells(j);
            TroubCells(j)=TroubCells(i);
            TroubCells(i)=amin;
        end
    end
end
TroubCells=TroubCells(1:Troubn);
