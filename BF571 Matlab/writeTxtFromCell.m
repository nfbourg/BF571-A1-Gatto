function writeTxtFromCell(cell,filename)

if nargin < 2
    fid = fopen('myfile.txt','w');
    disp('File savet at myfile.txt');
else
    fid = fopen(filename,'w');
end

for i = 1:length(cell)
    fprintf(fid,'%s\n',cell{i});
end
fclose(fid);