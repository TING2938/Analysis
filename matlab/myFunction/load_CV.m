function data = load_CV(fileName)
fid = fopen(fileName);
ii = 0;
while (~feof(fid))
   line = fgetl(fid);
   if strcmp(strtrim(line), 'End Comments')
       break;
   end
   ii = ii + 1;
   if ii > 100
       fseek(fid, 0, 'bof');
       break;
   end
end

data = cell2mat(textscan(fid, "%f %f %f%*[^\n]", "CommentStyle", "#"));
end