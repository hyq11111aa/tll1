function [Mesh]=Coarsenx(Mesh,Coarsnum, CoarsCs)
Cells=Mesh.Cells;   Cells_num=Mesh.Cellsn;
Ex=Mesh.Ex; Ex_num=Mesh.Exn;   Ey=Mesh.Ey; Ey_num=Mesh.Eyn;

k=1;  Ps=0;  Cellsp=zeros(1,Coarsnum);
while k<=Coarsnum
    Cp=CoarsCs(k);  Cpn=Cells(Cp).Parent;    Cp=Cells(Cpn);   cls=Cp.Cs;
    C1=Cells(cls(1));   C2=Cells(cls(2));    C3=Cells(cls(3));   C4=Cells(cls(4));
    Cp.U(:,1)=1/4*(C1.U(:,1)+C2.U(:,1)+C3.U(:,1)+C4.U(:,1));
    Cp.U(:,2)=3/8*(-C1.U(:,1)+C2.U(:,1)+C3.U(:,1)-C4.U(:,1))+...
        1/8*(C1.U(:,2)+C2.U(:,2)+C3.U(:,2)+C4.U(:,2));
    Cp.U(:,3)=3/8*(-C1.U(:,1)-C2.U(:,1)+C3.U(:,1)+C4.U(:,1))+...
        1/8*(C1.U(:,3)+C2.U(:,3)+C3.U(:,3)+C4.U(:,3));
    Ps=Ps+1;    Cellsp(Ps)=Cpn;     Cells(Cpn)=Cp;
    k=k+4;
end

ExN=Ex_num;  EyN=Ey_num;  CellsN=Cells_num;
Exn=1:Ex_num;  Eyn=1:Ey_num;  Cellsn=1:Cells_num;
for l=1:Ex_num
    Ex(l).ref=l;
end
for l=1:Ey_num
    Ey(l).ref=l;
end
for l=1:Cells_num
    Cells(l).ref=l;
end
    
for k=1:Ps
    Cpn=Cellsp(k);   Cc=Cells(Cpn);   cls=Cc.Cs;  eds=Cc.Edges;
    Cls1=Cells(cls(1));     Cls3=Cells(cls(3));
    EL=Ey(eds(1));  ER=Ey(eds(2));  EB=Ex(eds(3));  ET=Ex(eds(4)); 
    CL=Cells(EL.K1);    CR=Cells(ER.K2);    CB=Cells(EB.K1);    CT=Cells(ET.K2);

    if  CB.HaveChild 
        %%%%%%%%%%%%%%%%%%% update relation of edge and cells %%%%%%%%%
        CB4=Cells(CB.Cs(4)); Ex(CB4.Edges(4)).K2=Cpn;
        CB3=Cells(CB.Cs(3)); Ex(CB3.Edges(4)).K2=Cpn;
    else
        %%%%%%%%%%%%%%%%%%%% delete small edges in necessary %%%
        Ex(Cc.Edges(3)).HaveChild=0;    Ex(Cc.Edges(3)).Cs=zeros(1,2);
        ee=Cls1.Edges(3);   ee_newN=Ex(ee).ref;     
        Exn(ee_newN:ExN-2)=Exn(ee_newN+2:ExN);  ExN=ExN-2;
        for l=ee+2:Ex_num
            Ex(l).ref=Ex(l).ref-2;
        end
    end

    if  CL.HaveChild
        %%%%%%%%%%%%%%%%%%% update relation of edge and cells %%%%%%%%%
        CL2=Cells(CL.Cs(2)); Ey(CL2.Edges(2)).K2=Cpn;
        CL3=Cells(CL.Cs(3)); Ey(CL3.Edges(2)).K2=Cpn;
    else
        %%%%%%%%%%%%%%%%%%%% delete small edges in necessary %%%%%
        Ey(Cc.Edges(1)).HaveChild=0;    Ey(Cc.Edges(1)).Cs=zeros(1,2);
        ee=Cls1.Edges(1);  ee_newN=Ey(ee).ref;     
        Eyn(ee_newN:EyN-2)=Eyn(ee_newN+2:EyN);  EyN=EyN-2;
        for l=ee+2:Ey_num
            Ey(l).ref=Ey(l).ref-2;
        end 
    end
    
    %%%%%%%%%%%%%%%%%%%% delete nodes, small edges in necessary %%%%%
    ee=Cls1.Edges(2);   ee_newN=Ey(ee).ref;     
    Eyn(ee_newN:EyN-2)=Eyn(ee_newN+2:EyN);   EyN=EyN-2;
    for l=ee+2:Ey_num
        Ey(l).ref=Ey(l).ref-2;
    end
    ee=Cls1.Edges(4);   ee_newN=Ex(ee).ref;    
    Exn(ee_newN:ExN-2)=Exn(ee_newN+2:ExN); 	 ExN=ExN-2;
    for l=ee+2:Ex_num
        Ex(l).ref=Ex(l).ref-2;
    end
    
    if CR.HaveChild
        %%%%%%%%%%%%%%%%%%% update relation of edge and cells %%%%%%%%%
        CR1=Cells(CR.Cs(1)); Ey(CR1.Edges(1)).K1=Cpn;
        CR4=Cells(CR.Cs(4)); Ey(CR4.Edges(1)).K1=Cpn;
    else
        %%%%%%%%%%%%%%%%%%%% delete  small edges in necessary %%%%%
        Ey(Cc.Edges(2)).HaveChild=0; Ey(Cc.Edges(2)).Cs=zeros(1,2);
        ee=Cls3.Edges(2)-1;   ee_newN=Ey(ee).ref;    
        Eyn(ee_newN:EyN-2)=Eyn(ee_newN+2:EyN);  EyN=EyN-2;
        for l=ee+2:Ey_num
            Ey(l).ref=Ey(l).ref-2;
        end 
    end
    
    if CT.HaveChild
        %%%%%%%%%%%%%%%%%%% update relation of edge and cells %%%%%%%%%
        CT1=Cells(CT.Cs(1));    Ex(CT1.Edges(3)).K1=Cpn;
        CT2=Cells(CT.Cs(2));    Ex(CT2.Edges(3)).K1=Cpn;
    else
        %%%%%%%%%%%%%%%%%%%% delete small edges in necessary %%%%%
        Ex(Cc.Edges(4)).HaveChild=0; Ex(Cc.Edges(4)).Cs=zeros(1,2);
        ee=Cls3.Edges(4)-1;  ee_newN=Ex(ee).ref;     
        Exn(ee_newN:ExN-2)=Exn(ee_newN+2:ExN);   ExN=ExN-2;
        for l=ee+2:Ex_num
            Ex(l).ref=Ex(l).ref-2;
        end 
    end
    
    Cc=Cells(Cpn).Cs(1);         Cc_newN=Cells(Cc).ref; 
    Cells(Cpn).HaveChild=0;      Cells(Cpn).Cs=zeros(1,4);
    Cellsn(Cc_newN:CellsN-4)=Cellsn(Cc_newN+4:CellsN);  CellsN=CellsN-4;    
    for l=Cc+4:Cells_num
        Cells(l).ref=Cells(l).ref-4;
    end
end

for k=1:ExN
    Exc=Ex(Exn(k));  Exc.K1=Cells(Exc.K1).ref;  Exc.K2=Cells(Exc.K2).ref; 
    if Exc.HaveChild
        Exc.Cs(1)=Ex(Exc.Cs(1)).ref;   Exc.Cs(2)=Ex(Exc.Cs(2)).ref; 
    end
    if Exc.Parent>0
        Exc.Parent=Ex(Exc.Parent).ref;     
    end
    Ex(Exn(k))=Exc;
end
for k=1:EyN
    Eyc=Ey(Eyn(k));  Eyc.K1=Cells(Eyc.K1).ref;  Eyc.K2=Cells(Eyc.K2).ref;
    if Eyc.HaveChild
        Eyc.Cs(1)=Ey(Eyc.Cs(1)).ref;   Eyc.Cs(2)=Ey(Eyc.Cs(2)).ref; 
    end
    if Eyc.Parent>0
        Eyc.Parent=Ey(Eyc.Parent).ref;     
    end
    Ey(Eyn(k))=Eyc;
end
eds=zeros(1,4);
for k=1:CellsN
    Cc=Cells(Cellsn(k));    
    eds(1)=Ey(Cc.Edges(1)).ref;    eds(2)=Ey(Cc.Edges(2)).ref;
    eds(3)=Ex(Cc.Edges(3)).ref;    eds(4)=Ex(Cc.Edges(4)).ref;
    Cc.Edges=eds; 
    if Cc.HaveChild
        Cc.Cs(1)=Cells(Cc.Cs(1)).ref;   Cc.Cs(2)=Cells(Cc.Cs(2)).ref;
        Cc.Cs(3)=Cells(Cc.Cs(3)).ref;   Cc.Cs(4)=Cells(Cc.Cs(4)).ref;
    end
    if Cc.Parent>0
        Cc.Parent=Cells(Cc.Parent).ref;
    end
    Cells(Cellsn(k))=Cc;
end

Ex=Ex(Exn(1:ExN));   Ey=Ey(Eyn(1:EyN));   Cells=Cells(Cellsn(1:CellsN));
Ex_num=ExN;     Ey_num=EyN;     Cells_num=CellsN;
for l=1:Ex_num
    Ex(l).Indx=l;
end
for l=1:Ey_num
    Ey(l).Indx=l;
end
for l=1:Cells_num
    Cells(l).Indx=l;    Cells(l).NoTroubleC=1;
end

% ExN=Ex_num;    EyN=Ey_num;    CellsN=Cells_num;
% Exn=1:Ex_num;  Eyn=1:Ey_num;  Cellsn=1:Cells_num;
% for l=1:Ex_num
%     Ex(l).ref=l;
% end
% for l=1:Ey_num
%     Ey(l).ref=l;
% end
% for l=1:Cells_num
%     Cells(l).ref=l;
% end
% 
% for k=1:Parentn
%     Cpn=Cellsp(k); Cp=Cells(Cpn);   cls=Cp.Cs;  eds=Cp.Edges;
%     Cls1=Cells(cls(1));     Cls3=Cells(cls(3));
%     %%%%%%%%%%%%%%%%%%%%%%%%% disp('delete bottom edges')
%     EB=Ex(eds(3));  CB=Cells(EB.K1);
%     if CB.HaveChild
%         %%%%%%%%%%%%%%%%%%% update relation of edge and cells %%%%%%%%%
%         CB4=Cells(CB.Cs(4)); Ex(CB4.Edges(4)).K2=Cpn;
%         CB3=Cells(CB.Cs(3)); Ex(CB3.Edges(4)).K2=Cpn;
%     else
%       %%%%%%%%%%%%%%%%%%%% delete small edges in necessary %%%
%        Ex(eds(3)).HaveChild=0;    Ex(eds(3)).Cs=[0,0];
%        ee=Cls1.Edges(3);   ee_newN=Ex(ee).ref;     
%        Exn(ee_newN:ee_newN+1)=[];  ExN=ExN-2; 
%        for l=ee+2:Ex_num
%            Ex(l).ref=Ex(l).ref-2;
%        end
%     end
%   
%     EL=Ey(eds(1));  CL=Cells(EL.K1);
%     if CL.HaveChild
%         CL2=Cells(CL.Cs(2)); Ey(CL2.Edges(2)).K2=Cpn;
%         CL3=Cells(CL.Cs(3)); Ey(CL3.Edges(2)).K2=Cpn;
%     else 
%         Ey(eds(1)).HaveChild=0;    Ey(eds(1)).Cs=zeros(1,2);
%         ee=Cls1.Edges(1);  ee_newN=Ey(ee).ref;     
%         Eyn(ee_newN:ee_newN+1)=[];  EyN=EyN-2;
%         for l=ee+2:Ey_num
%             Ey(l).ref=Ey(l).ref-2;
%         end 
%     end
%     
%     %%%%%%%%%%%%%%%%%%%% delete small edges inside the element %%%%%
%     ee=Cls1.Edges(2);   ee_newN=Ey(ee).ref;     
%     Eyn(ee_newN:ee_newN+1)=[];  EyN=EyN-2;   
%     for l=ee+2:Ey_num
%         Ey(l).ref=Ey(l).ref-2;
%     end
% 
%     ee=Cls1.Edges(4);   ee_newN=Ex(ee).ref;    
%     Exn(ee_newN:ee_newN+1)=[]; 	ExN=ExN-2;  
%     for l=ee+2:Ex_num
%         Ex(l).ref=Ex(l).ref-2;
%     end
%     
%     %%%%%%%%%%%%%%%%%%%% delete small edges inside the element %%%%%
%     ER=Ey(eds(2));  CR=Cells(ER.K2);
%     if CR.HaveChild
%         CR1=Cells(CR.Cs(1)); Ey(CR1.Edges(1)).K1=Cpn;
%         CR4=Cells(CR.Cs(4)); Ey(CR4.Edges(1)).K1=Cpn;
%     else
%         Ey(eds(2)).HaveChild=0;    Ey(eds(2)).Cs=zeros(1,2);
%         ee=Cls3.Edges(2)-1;  ee_newN=Ey(ee).ref;     
%         Eyn(ee_newN:ee_newN+1)=[];  EyN=EyN-2; 
%         for l=ee+2:Ey_num
%             Ey(l).ref=Ey(l).ref-2;
%         end 
%     end
%    
%   %%%%%%%%%%%%%%%%%%%% delete small edges in necessary %%%
%     ET=Ex(eds(4));  CT=Cells(ET.K2);
%     if CT.HaveChild
%         CT1=Cells(CT.Cs(1)); Ex(CT1.Edges(3)).K1=Cpn;
%         CT2=Cells(CT.Cs(2)); Ex(CT2.Edges(3)).K1=Cpn;
%     else
%        Ex(eds(4)).HaveChild=0;    Ex(eds(4)).Cs=[0,0];
%        ee=Cls3.Edges(4)-1;   ee_newN=Ex(ee).ref;     
%        Exn(ee_newN:ee_newN+1)=[];  ExN=ExN-2; 
%        for l=ee+2:Ex_num
%            Ex(l).ref=Ex(l).ref-2;
%        end
%     end
%     
%     Cp=Cells(Cpn).Cs(1);         Cc_newN=Cells(Cp).ref; 
%     Cells(Cpn).HaveChild=0;      Cells(Cpn).Cs=zeros(1,4);
%     Cellsn(Cc_newN:CellsN-4)=Cellsn(Cc_newN+4:CellsN);  CellsN=CellsN-4;   
%     for l=Cp+4:Cells_num
%         Cells(l).ref=Cells(l).ref-4;
%     end
% end
% 
% for k=1:ExN
%     Exc=Ex(Exn(k));  ek1=Exc.K1;    ek2=Exc.K2;
%     Exc.K1=Cells(ek1).ref;       Exc.K2=Cells(ek2).ref;
%     if Exc.HaveChild
%         Exc.Cs(1)=Ex(Exc.Cs(1)).ref;   Exc.Cs(2)=Ex(Exc.Cs(2)).ref; 
%     end
%     if Exc.Parent>0
%         Exc.Parent=Ex(Exc.Parent).ref;  
%     end
%     Ex(Exn(k))=Exc;
% end
% 
% for k=1:EyN
%     Eyc=Ey(Eyn(k));   ek1=Eyc.K1;    ek2=Eyc.K2;
%     Eyc.K1=Cells(ek1).ref;     Eyc.K2=Cells(ek2).ref;
%     if Eyc.HaveChild
%         c1=Eyc.Cs(1);   c2=Eyc.Cs(2);
%         Eyc.Cs(1)=Ey(c1).ref;   Eyc.Cs(2)=Ey(c2).ref; 
%     end
%     cp=Eyc.Parent;
%     if cp>0
%         Eyc.Parent=Ey(cp).ref;     
%     end
%     Ey(Eyn(k))=Eyc;
% end
% 
% for k=1:CellsN
%     Cp=Cells(Cellsn(k));   edss=Cp.Edges;   cls=Cp.Cs;
%     Cp.Edges=[Ey(edss(1)).ref, Ey(edss(2)).ref, Ex(edss(3)).ref, Ex(edss(4)).ref]; 
%     if Cp.HaveChild
%         cls(1)=Cells(cls(1)).ref;   cls(2)=Cells(cls(2)).ref;
%         cls(3)=Cells(cls(3)).ref;   cls(4)=Cells(cls(4)).ref;
%     end
%     Cp.Cs=cls;
%     if Cp.Parent>0
%         Cp.Parent=Cells(Cp.Parent).ref;
%     end
%     Cells(Cellsn(k))=Cp;
% end
% 
% Ex=Ex(Exn(1:ExN));   Ey=Ey(Eyn(1:EyN));   Cells=Cells(Cellsn(1:CellsN));
% Ex_num=ExN;     Ey_num=EyN;     Cells_num=CellsN; 
% for l=1:Ex_num
%     Ex(l).Indx=l;
% end
% for l=1:Ey_num
%     Ey(l).Indx=l;
% end
% for l=1:Cells_num
%     Cells(l).Indx=l;   Cells(l).NoTroubleC=1;
% end
Mesh.Cells=Cells;   Mesh.Cellsn=Cells_num;
Mesh.Ex=Ex; Mesh.Exn=Ex_num;   Mesh.Ey=Ey; Mesh.Eyn=Ey_num;