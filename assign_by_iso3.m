
function T = assign_by_iso3(iso3, vals, T, varname)

if iscell(iso3)
    iso3 = string(iso3);
end

nrows = height(T);
vArr = zeros(nrows,1);

for i = 1:nrows
    Tiso3i = T.iso3{i};
    % iso3idx = find(iso3 == Tiso3i);
    vArr(i) = vals(iso3 == Tiso3i);
end

T.(varname) = vArr;

end