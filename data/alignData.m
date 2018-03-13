%% Load
load ./allTS.mat

%% Align
for i=1:size(allEMG,1)
    for j=1:size(allEMG,2)
        alignedEMG{i,j}=allEMG{i,j}.align(allEvents{i,j},{'FHS','STO','SHS','FTO'},16*[1,2,1,2]);
    end
end

%% Save
save alignedEMG.mat alignedEMG