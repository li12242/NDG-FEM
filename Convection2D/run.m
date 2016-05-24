ele = [20, 40, 60, 80, 100];
deg = 1:3;

for ie = 1:numel(ele)
    for id = 1:numel(deg)
        Convection2DSetUp(deg(id), ele(ie));
    end% for
end% for