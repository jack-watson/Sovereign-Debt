function fullnames = iso32fullname(iso2nameTbl, iso3)

% Takes in a list of iso3 codes (string array) and a table mapping iso3
% codes to country names; outputs list of corresponding country names
% (string array).

if iscell(iso3)
    iso3 = string(iso3);
end

ncodes = length(iso3);
fullnames = strings([ncodes,1]);

for i = 1:ncodes
    iso3i = iso3(i);
    tblIdx = find(iso2nameTbl.iso3 == iso3i);
    fullnames(i) = iso2nameTbl.name(tblIdx); 
end






end