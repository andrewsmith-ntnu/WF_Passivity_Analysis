function S = setNestedField(S, fieldName, value)
% setNestedField  Set nested struct field using dot-separated name
%
% S = setNestedField(S, 'a.b.c', value)
%
% This utility uses subsasgn to assign a nested field specified as a
% dot-separated string. It creates the appropriate substruct indexing.

parts = strsplit(fieldName, '.');
subs = struct('type', {}, 'subs', {});
for i = 1:length(parts)
    subs(i).type = '.';
    subs(i).subs = parts{i};
end
S = subsasgn(S, subs, value);
end
