function x = read_numeric_file(path)
%READ_NUMERIC_FILE Load a column of numbers from a text file.
% Supports European decimal comma.
txt = fileread(path);
% Normalize decimal comma to dot
txt = regexprep(txt, ',', '.');
% Split by whitespace/newlines
C = regexp(txt, '[\r\n\t ]+', 'split');
C = C(~cellfun('isempty',C));
x = str2double(C);
x = x(~isnan(x));
end

