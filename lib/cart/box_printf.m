% box_printf(string,box_char,box_len) prints a message sandwitched between dashed
% lines. The character length of the dashed lines is box_length. For
% example
%
%           box_printf('Hello','-',100)
%
% produces 
%
% ----------------------------------------------------------------------------------------------------
% Hello
% ----------------------------------------------------------------------------------------------------

function box_printf(message_str,box_char,box_len)


fprintf('\n');

for i = 1 : box_len
    
    fprintf(box_char);
    
end

fprintf('\n');

fprintf(message_str);

fprintf('\n');

for i = 1 : box_len
    
    fprintf(box_char);
    
end

fprintf('\n');

