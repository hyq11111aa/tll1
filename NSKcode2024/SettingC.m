function [Mesh, Mesh_fine]=SettingC(Mesh)
global mesh_fine
Cells=Mesh.Cells;   Cellsn=Mesh.Cellsn;
Eyn=Mesh.Eyn;    Ey=Mesh.Ey;    Exn=Mesh.Exn;    Ex=Mesh.Ex;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% finest cells %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fCells=zeros(1,Cellsn);  fCn=0; 
for i=1:Cellsn
    if Cells(i).HaveChild
        continue;
    end
    fCn=fCn+1;    Cells(i).Indx_fine=fCn;    fCells(fCn)=i;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% finest interior edges %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
FEx=zeros(1,Exn);   FEy=zeros(1,Eyn);
fExn=0;         fEyn=0;  
for i=1:Exn
    if Ex(i).HaveChild
        continue;
    end
    fExn=fExn+1;    Ex(i).Indx_fine=fExn;   FEx(fExn)=i;
end
for i=1:Eyn
    if Ey(i).HaveChild
        continue;
    end
    fEyn=fEyn+1;    Ey(i).Indx_fine=fEyn;   FEy(fEyn)=i;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% look for neighbor %%%%%%%%%%%%%%%%%%%%%
for i=1:fCn
    Cc=Cells(fCells(i));   Neighbn=0;  Neighb=zeros(1,8);  eds=Cc.Edges;
        
    CLn=Ey(eds(1)).K1; CL=Cells(CLn);
    if CL.HaveChild
        Neighb(Neighbn+1)=CL.Cs(2);   Neighb(Neighbn+2)=CL.Cs(3);  
        Neighbn= Neighbn+2;
    else
        Neighb(Neighbn+1)=CLn;    Neighbn= Neighbn+1;
    end
      
    CRn=Ey(eds(2)).K2;  CR=Cells(CRn);
    if CR.HaveChild
        Neighb(Neighbn+1)=CR.Cs(1);   Neighb(Neighbn+2)=CR.Cs(4);  
        Neighbn= Neighbn+2;
    else
        Neighb(Neighbn+1)=CRn;  Neighbn= Neighbn+1;
    end
    
    CBn=Ex(eds(3)).K1;    CB=Cells(CBn);
    if CB.HaveChild
        Neighb(Neighbn+1)=CB.Cs(4);   Neighb(Neighbn+2)=CB.Cs(3);  
        Neighbn= Neighbn+2;
    else
        Neighb(Neighbn+1)=CBn;  Neighbn= Neighbn+1;
    end

    CTn=Ex(eds(4)).K2;   CT=Cells(CTn);
    if CT.HaveChild
        Neighb(Neighbn+1)=CT.Cs(1);   Neighb(Neighbn+2)=CT.Cs(2);  
        Neighbn= Neighbn+2;
    else
        Neighb(Neighbn+1)=CTn;  Neighbn= Neighbn+1;
    end

    Cc.Neighb_num=Neighbn;  Cc.Neighb=Neighb;   Cells(fCells(i))=Cc;
end

Mesh.Cells=Cells;   Mesh.Ex=Ex; Mesh.Ey=Ey;

Mesh_fine=mesh_fine;
Mesh_fine.fEx=FEx;  Mesh_fine.fExn=fExn;   
Mesh_fine.fEy=FEy;  Mesh_fine.fEyn=fEyn;
Mesh_fine.fCells=fCells; Mesh_fine.fCellsn=fCn;

% figure
% hold on
% for k=1:fCn
%     Ck=Cells(fCells(k));    Levk=Ck.Lev;   center=Ck.Center;  h=Ck.Wid;   
%     xo=center(1);   yo=center(2);
%     xs=[xo-h/2, xo+h/2, xo+h/2, xo-h/2, xo-h/2];
%     ys=[yo-h/2, yo-h/2, yo+h/2, yo+h/2, yo-h/2];
%     if Levk==1
%         fill(xs, ys, 'b')
%     elseif Levk==2
%         fill(xs, ys, 'r')
%     elseif Levk==3
%         fill(xs, ys, 'g')
%     elseif Levk==4
%         fill(xs, ys, 'k')
%     else
%         fill(xs, ys, 'w')
%     end
% end
% 
