fid1=fopen('D:\tjuhanlu\NDG-FEM\NDG-FEM\SWE2D\@swe2d_obliquejump\mesh\quad1900\ObliqueHydraulicJump_Quad.msh');
fid2=fopen('D:\tjuhanlu\exercise of matlab\mesh.node','w');
fid3=fopen('D:\tjuhanlu\exercise of matlab\mesh.edge','w');
fid4=fopen('D:\tjuhanlu\exercise of matlab\mesh.ele','w');
fid5=fopen('D:\tjuhanlu\exercise of matlab\judge','w');
for i=1:4
   fgetl(fid1);
end
A=fscanf(fid1,'%s',1);
fprintf(fid2,'%s %d',A,2);
for i=1:str2num(A)+1
      temp=fgetl(fid1);
      fprintf(fid2,'%s\n',temp);
end

for i=1:2
   fgetl(fid1);
end
B=fscanf(fid1,'%d',1);
edge=0;
for i=1:B+1
      [temp]=fscanf(fid1,'%d %d %d %d %d %d %d\n',[1,7]);
      fprintf(fid5,'%d\n',temp(2));
      if temp(2)==1
          edge=edge+1;
      else
          m=temp(2);
          break;
      end
end
fseek(fid1,0,-1);
fprintf(fid3,'%d %d %d\n',edge,2,1);
C=str2num(A)+8;
for i=1:C
      fgetl(fid1);
end
for i=1:edge
    [temp]=fscanf(fid1,'%d %d %d %d %d %d %d\n',[1,7]);
    fprintf(fid3,'%d %d %d %d\n',temp(1),temp(6),temp(7),temp(4));
end

ele=B-edge;
fprintf(fid4,'%d %d %d\n',ele,4,3);
for i=1:ele
    if m == 3
        [temp]=fscanf(fid1,'%d %d %d %d %d %d %d %d %d\n',[1,9]);
        fprintf(fid4,'%d %d %d %d %d %d\n',i,temp(6),temp(7),temp(9),temp(8),temp(4));
    else
        [temp]=fscanf(fid1,'%d %d %d %d %d %d %d %d \n',[1,8]);
        fprintf(fid4,'%d %d %d %d %d\n',i,temp(6),temp(7),temp(8),temp(4));
    end
end

