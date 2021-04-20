function data = loadtxt(fileName, comments)
%% function loadtxt
%  read matrix-like ascii file, lines with comment style will be ignored.
if length(fileName) > 4 && strcmp(fileName(end-3:end), '.mat')
    data = load(fileName);
    return;
end

if (nargin < 2)
    comments = '#@'; % default comment style;
end
fid = fopen(fileName);
count = 0;
while (~feof(fid))
    line = fgetl(fid);
    if isempty(line) || contains(comments,line(1))
        count = count + 1;
    else
        break;
    end
end
fclose(fid);

fid = fopen(fileName);
data = cell2mat(textscan(fid, '' , 'headerlines', count));
fclose(fid);

end

