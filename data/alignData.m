%% Load
load ./controlData.mat

%% Align
for i=1:size(allEMG,1)
    for j=1:size(allEMG,2)
        alignedEMGControls{i,j}=allEMG{i,j}.align(allEvents{i,j},{'RHS','LTO','LHS','RTO'},64*[1,2,1,2]);
    end
end
clear allEMG allEvents allTS

%% Now patients:
load ./patientData.mat

%% Align
for i=1:size(allEMG,1)
    for j=1:size(allEMG,2)
        alignedEMGPatients{i,j}=allEMG{i,j}.align(allEvents{i,j},{'RHS','LTO','LHS','RTO'},64*[1,2,1,2]);
    end
end

%% Save
save alignedEMGPatients.mat alignedEMGPatients 
save alignedEMGControls.mat alignedEMGControls
