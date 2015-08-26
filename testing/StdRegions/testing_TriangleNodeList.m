function testing_TriangleNodeList
nOrder = 3;
vOrder = [2,3,1];
facelist = StdRegions.Triangle.reorderTriNodeList(nOrder, vOrder);
result = [4,7,9,10,3,6,8,2,5,1]';
if any(facelist - result)
    error('reorder face list error')
end

vOrder = [2,1,3];
facelist = StdRegions.Triangle.reorderTriNodeList(nOrder, vOrder);
result = [4,3,2,1,7,6,5,9,8,10]';
if any(facelist - result)
    error('reorder face list error')
end
end