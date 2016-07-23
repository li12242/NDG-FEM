ele = 20;
deg = 2;

for ie = 1:numel(ele)
    for id = 1:numel(deg)
        Convection2DSetUp(deg(id), ele(ie));
    end% for
end% for