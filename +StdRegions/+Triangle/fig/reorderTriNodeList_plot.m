temp{1} = find(abs(s+1)<10e-10);
temp{2} = find(abs(r+s)<10e-10);
temp{3} = find(abs(r+1)<10e-10);
%%
figure; 
subplot(2,2,1); hold on; axis equal
for i = 1:3
    plot(r(temp{i}),s(temp{i}), '-')
end

plot(r, s, 'o')
for i = 1:numel(r)
    text(r(i)+0.05, s(i)+0.05, num2str(i))
end
VerStr = {'A', 'B', 'C'};
i = 1;
text(r(i)-0.1, s(i)-0.1, VerStr{1})
i = 4;
text(r(i)-0.1, s(i)-0.1, VerStr{2})
i = 10;
text(r(i)-0.1, s(i)-0.1, VerStr{3})
xlabel('next Element: local indicator')
set(gca,'YLim', [-1.4, 1.5], 'XLim', [-1.4, 1.5])

%%
subplot(2,2,2); hold on; axis equal
for i = 1:3
    plot(r(temp{i}),s(temp{i}), '-')
end

plot(r, s, 'o')
for i = 1:numel(r)
    text(r(i)+0.05, s(i)+0.05, num2str(i))
end
VerStr = {'C', 'B', 'A'};
i = 1;
text(r(i)-0.1, s(i)-0.1, VerStr{1})
i = 4;
text(r(i)-0.1, s(i)-0.1, VerStr{2})
i = 10;
text(r(i)-0.1, s(i)-0.1, VerStr{3})
xlabel('this Element: local indicator')
set(gca,'YLim', [-1.4, 1.5], 'XLim', [-1.4, 1.5])

%%
subplot(2,2,3); hold on; axis equal
for i = 1:3
    plot(r(temp{i}),s(temp{i}), '-')
end

plot(r, s, 'o')
for i = 1:numel(r)
    text(r(i)+0.05, s(i)+0.05, num2str(i))
end
VerStr = {'A', 'B', 'C'};
i = 1;
text(r(i)-0.1, s(i)-0.1, VerStr{1})
i = 4;
text(r(i)-0.1, s(i)-0.1, VerStr{2})
i = 10;
text(r(i)-0.1, s(i)-0.1, VerStr{3})
xlabel('next Element: local indicator')
set(gca,'YLim', [-1.4, 1.5], 'XLim', [-1.4, 1.5])

%%
subplot(2,2,4); hold on; axis equal
for i = 1:3
    plot(r(temp{i}),s(temp{i}), '-')
end

plot(r, s, 'o')
for i = 1:numel(r)
    text(r(i)+0.05, s(i)+0.05, num2str(facelist(i)))
end
VerStr = {'C', 'B', 'A'};
i = 1;
text(r(i)-0.1, s(i)-0.1, VerStr{1})
i = 4;
text(r(i)-0.1, s(i)-0.1, VerStr{2})
i = 10;
text(r(i)-0.1, s(i)-0.1, VerStr{3})
xlabel('this Element: next Element local indicator')
set(gca,'YLim', [-1.4, 1.5], 'XLim', [-1.4, 1.5])

