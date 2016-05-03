nele = [50, 100, 200, 400, 600, 800, 1000];
hd = [1e-8, 1e-4, 1e-2, 1];

for ih = 1:numel(hd)
    for ie = 1:numel(nele)
        [h ,q] = SWESetUp(nele(ie), hd(ih));
    end
end

hd = 1e-16;

for ih = 1:numel(hd)
    for ie = 1:numel(nele)
        [h ,q] = SWESetUp(nele(ie), hd(ih));
    end
end