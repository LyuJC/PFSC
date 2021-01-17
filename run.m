f=load('Handwritten_numerals.mat'); 
data=f.data;
label=f.labels;


para1=[1e-4 1e-2 1 100 1e4];
para2=[1e-4 1e-2 1 100 1e4];

%normalization
for i=1:size(data,1)
dist = max(max(data{i})) - min(min(data{i}));
m01 = (data{i} - min(min(data{i})))/dist;
data{i} = 2 * m01 - 1;
end

%PFSC
for i=1:length(para1)
    for j=1:length(para2)
        [result,FV,YY,SV,Wv]=PFSC(data,label,para1(i),para2(j));
        dlmwrite('Handwritten_numerals.txt',[para1(i) para2(j) result(1,:) result(2,:) result(3,:) result(4,:) result(5,:) result(6,:) ],'-append','delimiter','\t','newline','pc');
    end
end

