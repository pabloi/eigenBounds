clearvars
clc

Nsub=16;
dataDir='/datadump/rawData/Exp0001/matData/HPF30/';
condNames={'TM base','Adapt','Post'}; %Retrieving late data for these three conditions

for i=1:Nsub
    %Load
    aux=num2str(i/100);
    load([dataDir 'P00' aux(3:end) '.mat'])
    for k=1:3
        %Get last trial in condition
        t=expData.metaData.getTrialsInCondition(condNames{k});
        t=t(end); %Last trial
        t1=expData.data{t}.gaitEvents.Time(end)-5; %Avoiding very last 5 sec
        t0=t1-40; %40 secs of data.
        %Get EMG
        allEMG{i,k}=expData.data{t}.procEMGData.split(t0,t1);
        %Get events
        allEvents{i,k}=expData.data{t}.gaitEvents.split(t0,t1);
        %Get kinematics
        allTS{i,k}=expData.data{t}.markerData.split(t0,t1);
    
    end
end

%% Save
save patientData allEMG allEvents allTS