function [Mesh]=Time2Newton(Mesh, Mesh_fine, t, dt)
global ko 
fCells=Mesh_fine.fCells;    fCellsn=Mesh_fine.fCellsn;
Cells=Mesh.Cells;   Cellsit=Cells; 
for k=1:fCellsn
    id=fCells(k);   Uk=Cellsit(id).U; 
    Uk(1,:)=Uk(1,:)+0.0001; Cellsit(id).U=Uk; 
end
[Cellsit]=ComputeZ(Cellsit, Mesh, Mesh_fine);  
nk=10*ko; 

ik=0;   kmax=50;   
AM=Jacobx(Cellsit, Mesh, Mesh_fine, dt); 
[Lm, Um]=ilu(AM);
while ik<=kmax
    Fu=Newton_Fun(Cellsit, Mesh, Mesh_fine, t, dt); 
%     s=-AM\Fu;       
    [s, flag1]=gmres(AM, -Fu, 50, 1.0e-8, 100, Lm, Um); 
    if flag1~=0
        disp(flag2)
    end
    for k=1:fCellsn
        idk=fCells(k);   Uc=Cellsit(idk).U+reshape(s((k-1)*nk+1:k*nk),ko,10)';
        Cellsit(idk).U=Uc;  
    end
    if norm(s,inf)<0.0001
       break;
    end
    ik=ik+1;  %disp([ik, norm(s, inf)])
end
%     ssssssss
disp(['it=' num2str(ik)]);
Mesh.Cells=Cellsit;