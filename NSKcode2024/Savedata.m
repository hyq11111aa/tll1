function Savedata(Cells, fCells, fCellsn,filenam)

 C0n=0;  C1n=0;  C2n=0;  C3n=0;
 C0=zeros(1,fCellsn);   C1=zeros(1,fCellsn);
 C2=zeros(1,fCellsn);   C3=zeros(1,fCellsn);
 for k=1:fCellsn
     Cc=Cells(fCells(k));
     if Cc.Lev==0
          C0n=C0n+1;  C0(C0n)=fCells(k);
      elseif Cc.Lev==1
          C1n=C1n+1;  C1(C1n)=fCells(k);
      elseif Cc.Lev==2
          C2n=C2n+1;  C2(C2n)=fCells(k);
     else
          C3n=C3n+1;  C3(C3n)=fCells(k);
      end   
 end
 
fid=fopen(filenam,'wt');
fprintf(fid, '%s\n', 'TITLE="solution of the isothermal NSK eqs"');
fprintf(fid, '%s\n', 'VARIABLES= "x", "y", "rho"');

% figure 
% hold on
if C0n>0
    fprintf(fid, '%s %d %s %d\n', 'ZONE T="LEV=0", ZONETYPE=FEQUADRILATERAL, DATAPACKING=POINT, N=',...
        C0n*4, ', E=', C0n);
    for k=1:C0n
        Cc=Cells(C0(k));    xo=Cc.Center(1);    yo=Cc.Center(2);    hh=0.5*Cc.Wid;
        r=Cc.U(1,:);
        fprintf(fid, '%8.6f %8.6f %10.6f\n', xo-hh, yo-hh, r(1)-r(2)-r(3));
        fprintf(fid, '%8.6f %8.6f %10.6f\n', xo+hh, yo-hh, r(1)+r(2)-r(3));
        fprintf(fid, '%8.6f %8.6f %10.6f\n', xo+hh, yo+hh, r(1)+r(2)+r(3));
        fprintf(fid, '%8.6f %8.6f %10.6f\n', xo-hh, yo+hh, r(1)-r(2)+r(3));
%         xl=xo-hh;   xr=xo+hh;   yb=yo-hh;   yt=yo+hh;
%         fill([xl xr xr xl xl], [yb yb yt yt yb],'w')
    end
    for k=1:C0n
        fprintf(fid, '%d %d %d %d\n', (k-1)*4+1,(k-1)*4+2, (k-1)*4+3, k*4);
    end
end

if C1n>0
    fprintf(fid, '%s %d %s %d\n', 'ZONE T="LEV=1", ZONETYPE=FEQUADRILATERAL, DATAPACKING=POINT, N=',...
        C1n*4, ', E=', C1n);
    for k=1:C1n
        Cc=Cells(C1(k));    xo=Cc.Center(1);    yo=Cc.Center(2);    hh=0.5*Cc.Wid;
        r=Cc.U(1,:);
        fprintf(fid, '%8.6f %8.6f %10.6f\n', xo-hh, yo-hh, r(1)-r(2)-r(3));
        fprintf(fid, '%8.6f %8.6f %10.6f\n', xo+hh, yo-hh, r(1)+r(2)-r(3));
        fprintf(fid, '%8.6f %8.6f %10.6f\n', xo+hh, yo+hh, r(1)+r(2)+r(3));
        fprintf(fid, '%8.6f %8.6f %10.6f\n', xo-hh, yo+hh, r(1)-r(2)+r(3));
%          xl=xo-hh;   xr=xo+hh;   yb=yo-hh;   yt=yo+hh;
%         fill([xl xr xr xl xl], [yb yb yt yt yb],'b')
    end
    for k=1:C1n
        fprintf(fid, '%d %d %d %d\n', (k-1)*4+1,(k-1)*4+2, (k-1)*4+3, k*4);
    end
end

if C2n>0
    fprintf(fid, '%s %d %s %d\n', 'ZONE T="LEV=2", ZONETYPE=FEQUADRILATERAL, DATAPACKING=POINT, N=',...
        C2n*4, ', E=', C2n);
    for k=1:C2n
        Cc=Cells(C2(k));    xo=Cc.Center(1);    yo=Cc.Center(2);    hh=0.5*Cc.Wid;
        r=Cc.U(1,:);
        fprintf(fid, '%8.6f %8.6f %10.6f\n', xo-hh, yo-hh, r(1)-r(2)-r(3));
        fprintf(fid, '%8.6f %8.6f %10.6f\n', xo+hh, yo-hh, r(1)+r(2)-r(3));
        fprintf(fid, '%8.6f %8.6f %10.6f\n', xo+hh, yo+hh, r(1)+r(2)+r(3));
        fprintf(fid, '%8.6f %8.6f %10.6f\n', xo-hh, yo+hh, r(1)-r(2)+r(3));
%         xl=xo-hh;   xr=xo+hh;   yb=yo-hh;   yt=yo+hh;
%         fill([xl xr xr xl xl], [yb yb yt yt yb],'r')
    end
    for k=1:C2n
        fprintf(fid, '%d %d %d %d\n', (k-1)*4+1,(k-1)*4+2, (k-1)*4+3, k*4);
    end
end

if C3n>0
    fprintf(fid, '%s %d %s %d\n', 'ZONE T="LEV=3", ZONETYPE=FEQUADRILATERAL, DATAPACKING=POINT, N=',...
        C3n*4, ', E=', C3n);

    for k=1:C3n
        Cc=Cells(C3(k));    xo=Cc.Center(1);    yo=Cc.Center(2);    hh=0.5*Cc.Wid;
        r=Cc.U(1,:);
        fprintf(fid, '%8.6f %8.6f %10.6f\n', xo-hh, yo-hh, r(1)-r(2)-r(3));
        fprintf(fid, '%8.6f %8.6f %10.6f\n', xo+hh, yo-hh, r(1)+r(2)-r(3));
        fprintf(fid, '%8.6f %8.6f %10.6f\n', xo+hh, yo+hh, r(1)+r(2)+r(3));
        fprintf(fid, '%8.6f %8.6f %10.6f\n', xo-hh, yo+hh, r(1)-r(2)+r(3));
%         xl=xo-hh;   xr=xo+hh;   yb=yo-hh;   yt=yo+hh;
%         fill([xl xr xr xl xl], [yb yb yt yt yb],'g')
    end
    for k=1:C3n
        fprintf(fid, '%d %d %d %d\n', (k-1)*4+1,(k-1)*4+2, (k-1)*4+3, k*4);
    end
end

fclose(fid);
disp('finish file')