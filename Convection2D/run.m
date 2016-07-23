ele      = 40;
deg      = 1;
meshtype = 'quad';

for ie = 1:numel(ele)
    for id = 1:numel(deg)
        Convection2DSetUp(meshtype, deg(id), ele(ie));
    end% for
end% for