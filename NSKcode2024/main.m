% A=load('B50cnL0t30.mat');   ME=A.ME;
% it=length(find(ME(:,1)>0));  T=30;
% B=load('B50cnL3t30.mat');   MEb=B.ME;
% itb=length(find(MEb(:,1)>0));
% figure
% hold on
% plot(ME(:,1),ME(:,2)-ME(1,2),'r-','LineWidth',2)
% plot(MEb(:,1),MEb(:,2)-MEb(1,2),'k-','LineWidth',2)
% set(gca,'fontsize',20);
% axis([0, T, -1.0e-12, 1.0e-12])
% xlabel('time','FontSize',20)
% ylabel('mass loss','FontSize',20)
% box on
% 
% figure
% hold on
% plot(ME(:,1),ME(:,3),'r-','LineWidth',2)
% plot(MEb(:,1),MEb(:,3),'k-','LineWidth',2)
% set(gca,'fontsize',20);
% xlabel('time','FontSize',20)
% ylabel('energy','FontSize',20)
% box on

% figure
% plot(ME(2:it-1,1),...
%     -(ME(3:it,3)-ME(1:it-2,3))./(ME(3:it,1)-ME(1:it-2,1)),'r-','LineWidth',2)
% hold on
% plot(MEb(2:itb-1,1),...
%     -(MEb(3:itb,3)-MEb(1:itb-2,3))./(MEb(3:itb,1)-MEb(1:itb-2,1)),'k-','LineWidth',2)
% set(gca,'fontsize',20);
% xlabel('time','FontSize',20)
% ylabel('-dEh/dt','FontSize',20)
% box on

% figure
% semilogy(ME(2:it-1,1),...
%     -(ME(3:it,3)-ME(1:it-2,3))./(ME(3:it,1)-ME(1:it-2,1)),'r-','LineWidth',2)
% hold on
% semilogy(MEb(2:itb-1,1),...
%     -(MEb(3:itb,3)-MEb(1:itb-2,3))./(MEb(3:itb,1)-MEb(1:itb-2,1)),'k-','LineWidth',2)
% set(gca,'fontsize',20);
% xlabel('time','FontSize',20)
% ylabel('-dEh/dt','FontSize',20)
% box on

% figure
% plot(ME(2:it-1,1),...
%     log10(-(ME(3:it,3)-ME(1:it-2,3))./(ME(3:it,1)-ME(1:it-2,1))),'r-','LineWidth',2)
% hold on
% plot(MEb(2:itb-1,1),...
%     log10(-(MEb(3:itb,3)-MEb(1:itb-2,3))./(MEb(3:itb,1)-MEb(1:itb-2,1))),'k-','LineWidth',2)
% set(gca,'fontsize',20);
% xlabel('time','FontSize',20)
% ylabel('-dEh/dt','FontSize',20)
% box on

% figure
% hold on
% % plot(ME(:,1),ME(:,4),'k-','LineWidth',2)
% plot(MEb(:,1),MEb(:,4),'k-','LineWidth',2)
% set(gca,'fontsize',20);
% xlabel('time','FontSize',20)
% ylabel('number of elements','FontSize',20)
% box on
% return

clear;  clc;  
MGset;
global cell edge  hc LEV
M=200;   T=30;   ko=3;  LEV=0;  hc=1/M;
eold1=1.0;    eold2=1.0;    eold3=1.0;  eold4=1.0;    eold5=1.0;    
cell=struct('Indx',0,'Indx_fine', 0, 'Lev',0, 'HaveChild',0, 'Wid', 0,...
    'NoTrouble', 1, 'Parent', 0, 'Edges', zeros(1,4),'Center',[0, 0],'Cs', zeros(1,4),...
    'Neighb_num', 0, 'Neighb', zeros(1,8),'U', zeros(10,ko),...
    'eta', 1,  'NoTroubleC', 1, 'Coarse',0, 'ref', 0);
edge=struct('Indx',0,'Indx_fine', 0, 'Lev', 0, 'HaveChild', 0, 'Parent', 0,...
    'K1', 0, 'K2', 0,'Cs', zeros(1,2),'ref', 0);
Es(1:M*M*50)=edge;  Eb(1:ceil(M*M/2))=edge; Cs(1:M*M*20)=cell;
mesh=struct('Cells', Cs, 'Cellsn',0,'Ex', Es, 'Exn', 0, 'Ey', Es, 'Eyn', 0);
idf=zeros(1, M*M*50); idbf=zeros(1, M*M);
mesh_fine=struct('fCells', idf,'fCellsn',0,'fEx', idbf, 'fExn', 0, 'fEy', idbf, 'fEyn', 0);

[Mesh, Mesh_fine]=Init(M, LEV, 0);
% G=load('B50cnL0t20.mat'); tstart=G.tc; itstart=G.it;
%  MEstart=G.ME; Mesh=G.Mesh; Mesh_fine=G.Mesh_fine;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%save data
% Cells=Mesh.Cells;  fCells=Mesh_fine.fCells;  fCellsn=Mesh_fine.fCellsn;
% filenam=['B50cnL' num2str(LEV) 't0.dat'];
% Savedata(Cells, fCells, fCellsn, filenam);
tstart=0;  itstart=0; MEstart=zeros(10^5,4);  
Tc=0; dTc=1;
while Tc<T
    if Tc+dTc>T
        dTc=T-Tc;
    end
    Tc=Tc+dTc;  
    [Mesh, Mesh_fine, tc, it, ME]=LDG_time2(Mesh, Mesh_fine, tstart, itstart, MEstart, Tc);
    tstart=tc;  itstart=it; MEstart=ME;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%save data
    Cells=Mesh.Cells;  fCells=Mesh_fine.fCells;  fCellsn=Mesh_fine.fCellsn;
    filenam=['B50cnL' num2str(LEV) 't' num2str(tc) '.dat'];
    Savedata(Cells, fCells, fCellsn, filenam);
    ME=ME(1:it,:);
    save(['B50cnL' num2str(LEV) 't' num2str(tc) '.mat'], 'tc', 'it', 'ME')
end

figure
plot(ME(:,1),ME(:,2)-ME(1,2),'k-','LineWidth',2)
set(gca,'fontsize',20);
axis([0, T, -1.0e-12, 1.0e-12])
xlabel('time','FontSize',20)
ylabel('mass loss','FontSize',20)

figure
plot(ME(:,1),ME(:,3),'k-','LineWidth',2)
set(gca,'fontsize',20);
xlabel('time','FontSize',20)
ylabel('energy','FontSize',20)

figure
semilogy(ME(2:it-1,1),-(ME(3:it,3)-ME(1:it-2,3))./(ME(3:it,1)-ME(1:it-2,1)),'k-','LineWidth',2)
set(gca,'fontsize',20);
xlabel('time','FontSize',20)
ylabel('-dEh/dt','FontSize',20)

figure
plot(ME(:,1),ME(:,4),'k-','LineWidth',2)
set(gca,'fontsize',20);
xlabel('time','FontSize',20)
ylabel('number of computing cells','FontSize',20)


sum(ME(:,4))/it
% data_unif=load('time15.mat'); 
% MEunif=data_unif.ME; figure;  hold on
% plot(MEunif(2:end,1),MEunif(2:end,3),'r-','LineWidth',2)
% set(gca,'fontsize',20);
% xlabel('time','FontSize',20)
% ylabel('energy','FontSize',20)
% 
% data_adap=load('time15ME.mat'); 
% MEadap=data_adap.MEv; 
% plot(MEadap(2:end,1),MEadap(2:end,3),'k--','LineWidth',2)
% set(gca,'fontsize',20);
% xlabel('time','FontSize',20)
% ylabel('energy','FontSize',20)