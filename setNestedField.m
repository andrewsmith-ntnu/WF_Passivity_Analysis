function S = setNestedField(S, fieldName, value)
    parts = strsplit(fieldName, '.');
    subs = struct('type', {}, 'subs', {});
    for i = 1:length(parts)
        subs(i).type = '.';
        subs(i).subs = parts{i};
    end
    S = subsasgn(S, subs, value);
end
