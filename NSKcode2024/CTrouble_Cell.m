function [CoarsCs, Coarsnum]=CTrouble_Cell(Mesh, Mesh_fine)

fCells=Mesh_fine.fCells;    fCellsn=Mesh_fine.fCellsn;  
Cells=Mesh.Cells;   Ey=Mesh.Ey;     Ex=Mesh.Ex;  
Ctroubn=0;  CoarseCells=zeros(1,floor(fCellsn/8)); eta2=0.005;
for k=1:fCellsn
    id=fCells(k);   Ck=Cells(id);
    if (Ck.eta < eta2) && (Ck.Parent > 0) 
        Ctroubn=Ctroubn+1;  CoarseCells(Ctroubn)=id;   Cells(id).NoTroubleC=0; 
    end
end
% %%%%%%%%%%%% sort%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:Ctroubn
    for j=i+1:Ctroubn
        if CoarseCells(j) > CoarseCells(i)
            amax=CoarseCells(j);  CoarseCells(j)=CoarseCells(i);   CoarseCells(i)=amax;
        end
    end
end

l=1; Coarsnum=0;     CoarsCs=zeros(1,Ctroubn);
while l<=Ctroubn
    ll=CoarseCells(l);  Cc=Cells(ll);   Cp=Cc.Parent;
    %%%%%%%%%%%%% all the children of Cp have coarsing flag %%%%%%%%%%%%
    Clds=Cells(Cp).Cs;   C1=Cells(Clds(1));  C2=Cells(Clds(2));  
    C3=Cells(Clds(3));   C4=Cells(Clds(4)); 
    if C1.NoTroubleC
        Cells(ll).NoTroubleC=1; l=l+1;  continue;
    end
    if C2.NoTroubleC
        Cells(ll).NoTroubleC=1; l=l+1;  continue;
    end
    if C3.NoTroubleC
        Cells(ll).NoTroubleC=1; l=l+1;  continue;
    end
    if C4.NoTroubleC
        Cells(ll).NoTroubleC=1; l=l+1;  continue;
    end
   %%%% assume level difference of Cp and its neighbor is at most 1 
    %%%%%%%%%%%%% neighbors of Cp %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%% C1 Left, Bottom %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    C1EL=Ey(C1.Edges(1));
    if C1EL.HaveChild
        Cells(ll).NoTroubleC=1; l=l+1;  continue;
    end
    C1EB=Ex(C1.Edges(3));
    if C1EB.HaveChild
        Cells(ll).NoTroubleC=1; l=l+1;  continue;
    end
    %%%%%%%%%%%%%%%%%%%%%%% C2 bottom, right %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    C2EB=Ex(C2.Edges(3));
    if C2EB.HaveChild
        Cells(ll).NoTroubleC=1; l=l+1;  continue; 
    end
    C2ER=Ey(C2.Edges(2));
    if C2ER.HaveChild
        Cells(ll).NoTroubleC=1; l=l+1;  continue;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%% C3 right, top %%%%%%%%%%%%%%%%%%%%%%%%%%
    C3ER=Ey(C3.Edges(2));
    if C3ER.HaveChild
        Cells(ll).NoTroubleC=1; l=l+1;  continue;
    end
    C3ET=Ex(C3.Edges(4));
    if C3ET.HaveChild
        Cells(ll).NoTroubleC=1; l=l+1;  continue;
    end
    %%%%%%%%%%%%%%%%%%%%%% C4 left, top %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    C4EL=Ey(C4.Edges(1));
    if C4EL.HaveChild
        Cells(ll).NoTroubleC=1; l=l+1;  continue;
    end
    C4ET=Ex(C4.Edges(4));
    if C4ET.HaveChild
       Cells(ll).NoTroubleC=1; l=l+1;  continue;
    end

    %%%%%%%%% for all neighbors Cn of Cp, Cn.Lev-Cp.Lev<=1, coarsen all
    %%%%%%%%% children of Cp to Cp
    CoarsCs(Coarsnum+1:Coarsnum+4)=ll-3:ll;
    Coarsnum=Coarsnum+4;        l=l+4;  
end
%%%%%%%%%%%% sort%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:Coarsnum
    for j=i+1:Coarsnum
        if CoarsCs(j) > CoarsCs(i)
            amax=CoarsCs(j);  CoarsCs(j)=CoarsCs(i);    CoarsCs(i)=amax;
        end
    end
end
CoarsCs=CoarsCs(1:Coarsnum);
% return
% %%%%%%%%%%%%%%%%%%%%%%% look for trouble cells %%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ctroubn=0;  CoarCells=zeros(1,floor(fCellsn/4));    eta2=0.005;
% for k=1:fCellsn 
%     id=fCells(k);   Ck=Cells(id); 
%     if Ck.Parent == 0
%         continue;
%     end
%     if (Ck.eta < eta2) && (Ck.eta >0) %%
%         Ctroubn=Ctroubn+1;  CoarCells(Ctroubn)=id;    Cells(id).NoTroubleC=0; 
%     end
% end
% %%%%%%%%%%%% sort%%%%%%%%%%%%%%%%%%%%%%%%%
% for i=1:Ctroubn
%     for j=i+1:Ctroubn
%         if CoarCells(j) > CoarCells(i)
%             amax=CoarCells(j);  CoarCells(j)=CoarCells(i);   CoarCells(i)=amax;
%         end
%     end
% end
% 
% l=1;     Coarsnum=0;     CoarsCs=zeros(1,Ctroubn);
% while l<=Ctroubn
%     ll=CoarCells(l);  Cc=Cells(ll);   Cp=Cells(Cc.Parent);  
%     %%%%%%%%%%%%% all the children of Cp have coarsing flag %%%%%%%%%%%%
%     Clds=Cp.Cs;   C1=Cells(Clds(1));  C2=Cells(Clds(2));  
%     C3=Cells(Clds(3));   C4=Cells(Clds(4)); 
%     if C1.NoTroubleC
%         Cells(ll).NoTroubleC=1; l=l+1;  continue;
%     end
%     if C2.NoTroubleC
%         Cells(ll).NoTroubleC=1; l=l+1;  continue;
%     end
%     if C3.NoTroubleC
%         Cells(ll).NoTroubleC=1; l=l+1;  continue;
%     end
%     if C4.NoTroubleC
%         Cells(ll).NoTroubleC=1; l=l+1;  continue;
%     end
%     
%    %%%%assume level difference of Cp and its neighbor is less than 2  
%     %%%%%%%%%%%%% neighbors of Cp %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %%%%%%%%%%%%%%%%%%%%%%%% C1 Left, Bottom %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%        %%%%%%%%%%%%%%%%%%%%%%%% C1 Left, Bottom %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     C1EL=Ey(C1.Edges(1));
%     if C1EL.HaveChild
%         Cells(ll).NoTroubleC=1; l=l+1;  continue;
%     end
%     C1EB=Ex(C1.Edges(3));
%     if C1EB.HaveChild
%         Cells(ll).NoTroubleC=1; l=l+1;  continue;
%     end
%     %%%%%%%%%%%%%%%%%%%%%%% C2 bottom, right %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     C2EB=Ex(C2.Edges(3));
%     if C2EB.HaveChild
%         Cells(ll).NoTroubleC=1; l=l+1;  continue; 
%     end
%     C2ER=Ey(C2.Edges(2));
%     if C2ER.HaveChild
%         Cells(ll).NoTroubleC=1; l=l+1;  continue;
%     end
%     %%%%%%%%%%%%%%%%%%%%%%%%% C3 right, top %%%%%%%%%%%%%%%%%%%%%%%%%%
%     C3ER=Ey(C3.Edges(2));
%     if C3ER.HaveChild
%         Cells(ll).NoTroubleC=1; l=l+1;  continue;
%     end
%     C3ET=Ex(C3.Edges(4));
%     if C3ET.HaveChild
%         Cells(ll).NoTroubleC=1; l=l+1;  continue;
%     end
%     %%%%%%%%%%%%%%%%%%%%%% C4 left, top %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     C4EL=Ey(C4.Edges(1));
%     if C4EL.HaveChild
%         Cells(ll).NoTroubleC=1; l=l+1;  continue;
%     end
%     C4ET=Ex(C4.Edges(4));
%     if C4ET.HaveChild
%        Cells(ll).NoTroubleC=1; l=l+1;  continue;
%     end
%     %%%%%%%%% for all neighbors Cn of Cp, Cn.Lev-Cp.Lev<=1, coarsen all
%     %%%%%%%%% children of Cp to Cp
%     CoarsCs(Coarsnum+1:Coarsnum+4)=ll-3:ll;
%     Coarsnum=Coarsnum+4;        l=l+4;  
% end
% 
% %%%%%%%%%%%% sort%%%%%%%%%%%%%%%%%%%%%%%%%
% for i=1:Coarsnum
%     for j=i+1:Coarsnum
%         if CoarsCs(j) > CoarsCs(i)
%             amax=CoarsCs(j);  CoarsCs(j)=CoarsCs(i);    CoarsCs(i)=amax;
%         end
%     end
% end
% CoarsCs=CoarsCs(1:Coarsnum);
% 
