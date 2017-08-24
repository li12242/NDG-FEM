function read_gmsh(msh_file_name)
msh_fp = fopen(msh_file_name);
[path_str, name, ext] = fileparts(msh_file_name);

node_fp=fopen( fullfile(path_str, [name,'.node']),'w' );
edge_fp=fopen( fullfile(path_str, [name,'.edge']),'w' );
ele_fp=fopen( fullfile(path_str, [name,'.ele']),'w' );
%fid5=fopen( fullfile(path_str, [name,'.judge']),'w' );

for i=1:4
   fgetl(msh_fp);
end
Nv = fscanf(msh_fp,'%d',1); % node number
fgetl(msh_fp);
fprintf(node_fp,'%d %d\n', Nv,2);
for i=1:(Nv)
      temp=fgetl(msh_fp);
      data = str2num(temp);
      fprintf(node_fp,'%d %g %g\n', data(1:3));
end

for i=1:2
   fgetl(msh_fp);
end
K = fscanf(msh_fp,'%d',1); % element number
edge=0;
for i=1:(K+1)
      [temp] = fscanf(msh_fp,'%d %d %d %d %d %d %d\n',[1,7]);
      %fprintf(fid5,'%d\n',temp(2));
      if temp(2) == 1
          edge=edge+1;
      else
          m = temp(2);
          break;
      end
end
fseek(msh_fp,0,-1);
fprintf(edge_fp,'%d %d %d\n',edge,2,1);
C=Nv+8;
for i=1:C
      fgetl(msh_fp);
end
for i=1:edge
    [temp]=fscanf(msh_fp,'%d %d %d %d %d %d %d\n',[1,7]);
    fprintf(edge_fp,'%d %d %d %d\n',temp(1),temp(6),temp(7),temp(4));
end

Ne=K-edge;
if m == 3
    fprintf(ele_fp,'%d %d %d\n',Ne,4,3);
else
    fprintf(ele_fp,'%d %d %d\n',Ne,3,3);
end
for i=1:Ne
    if m == 3
        [temp]=fscanf(msh_fp,'%d %d %d %d %d %d %d %d %d\n',[1,9]);
        fprintf(ele_fp,'%d %d %d %d %d %d\n',i,temp(6),temp(7),temp(9),temp(8),temp(4));
    else
        [temp]=fscanf(msh_fp,'%d %d %d %d %d %d %d %d \n',[1,8]);
        fprintf(ele_fp,'%d %d %d %d %d\n',i,temp(6),temp(7),temp(8),temp(4));
    end
end

fclose(msh_fp);
fclose(node_fp);
fclose(ele_fp);
fclose(edge_fp);
end% func
