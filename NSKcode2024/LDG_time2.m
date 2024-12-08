function  [Mesh, Mesh_fine, t, it, ME]=LDG_time2(Mesh, Mesh_fine, tstart, itstart, MEstart, Tc)
global hc LEV
t=tstart;   dt=0.1*hc/2^LEV;  it=itstart;  ME=MEstart;
MEs=MEfun(Mesh, Mesh_fine, t);  ME(it+1,:)=MEs;

while t<Tc
    if t+dt>=Tc
        dt=Tc-t;
    end
%     if rem(it+1,10)==0
%         [Mesh, TroubCells, Troubn]=Trouble_Cellx(Mesh, Mesh_fine, LEV);
%         Troubn
%         if Troubn>0
%             [Mesh]=Refine(Mesh, TroubCells, Troubn);
%             [Mesh, Mesh_fine]=SettingC(Mesh);
%         end
%         [CoarsCs, Coarsnum]=CTrouble_Cell(Mesh, Mesh_fine);
%         Coarsnum
%         if Coarsnum>0
%             [Mesh]=Coarsenx(Mesh, Coarsnum, CoarsCs);
%             [Mesh, Mesh_fine]=SettingC(Mesh);
%         end
%     end
    [Mesh]=Time2Newton(Mesh, Mesh_fine, t, dt);
        
    t=t+dt;  it=it+1; disp(['t=' num2str(t)] )
    
    MEs=MEfun(Mesh, Mesh_fine, t);  ME(it,:)=MEs;
end
ME=ME(1:it,:);