%% broken files
close all
clear all
clc

tasks2BeReNamed = {'20240307-133830-133839-GraspObject', ...
    '20240307-133830-134703-GraspObject', ...
    '20240307-135309-135446-GraspObject', ...
    '20240307-135309-140133-GraspObject', ...
    '20240307-140857-140906-GraspObject', ...
    '20240307-140857-141641-GraspObject'};

%% rename broken files
for ii = 1:length(tasks2BeReNamed)
    tasks2BeReNamed{ii};
    tempDataFile = load(['\\131.215.27.23\Raid01\Data\S4\20240307\Task\' tasks2BeReNamed{ii} '.mat']);
    hDebug = Debug.Debugger('loadDataFromFramework','screen');
    hArrayMPF = hst.Array('MPF', 'S4');
    hSubjectS4 = hst.Subject('S4', hDebug);
    NeuralSource.hArrayMaps = {hArrayMPF.hArrayMap};
    tempDataFile.NeuralSource.hArrayMaps = {hArrayMPF.hArrayMap};
    tempDataFile.Options.arrays          = {hArrayMPF};
    tempDataFile.Options.subject         = {hSubjectS4};
    tempDataFile.Options.output          = replace(tempDataFile.Options.output,'S3','S4'); 
    tempDataFile.NeuralSource.hXIPPMEX.arrayString = {hArrayMPF.ID};
    tempDataFile.NeuralSource.hXIPPMEX.arrayString = {hArrayMPF.ID};
    tempDataFile.NeuralSource.output = replace(tempDataFile.NeuralSource.output,'S3','S4');
    tempDataFile.NeuralSource.hCBMEX.outputPath = replace(tempDataFile.NeuralSource.hCBMEX.outputPath,'S3','S4');
    tempDataFile.NeuralSource.hXIPPMEX.recordFilenames =  replace(tempDataFile.NeuralSource.hXIPPMEX.recordFilenames,'MPx','MPF');
    tempDataFile.NeuralSource.hXIPPMEX.recordDirectories = replace(tempDataFile.NeuralSource.hXIPPMEX.recordDirectories,'MPx','MPF');
    tempDataFile.NeuralSource.hXIPPMEX.recordFilenames =  replace(tempDataFile.NeuralSource.hXIPPMEX.recordFilenames,'S3','S4');
    tempDataFile.NeuralSource.hXIPPMEX.recordDirectories = replace(tempDataFile.NeuralSource.hXIPPMEX.recordDirectories,'S3','S4');
    for tr = 1:Task.nTrials
        tempDataFile.Task.TrialData(tr).neu_filenames = replace(tempDataFile.Task.TrialData(tr).neu_filenames,'MPx','MPF');
    end
    % save
    idString        = tempDataFile.idString;
    Runtime         = tempDataFile.Runtime;
    Data            = tempDataFile.Data;
    Options         = tempDataFile.Options;
    NeuralSource    = tempDataFile.NeuralSource;
    Sync            = tempDataFile.Sync;
    Task            = tempDataFile.Task;
%     save(['\\131.215.27.23\Raid01\Data\S4\20240206\Task\' tasks2BeReNamed{ii}], 'idString', ...
%                                                                                 'Runtime', ...
%                                                                                 'Data', ...
%                                                                                 'Options', ...
%                                                                                 'NeuralSource');
    save(['\\131.215.27.23\Raid01\Data\S4\20240307\Task\' tasks2BeReNamed{ii}], 'idString', ...
                                                                                'Runtime', ...
                                                                                'Data', ...
                                                                                'Options', ...
                                                                                'NeuralSource', ...
                                                                                'Sync', ...
                                                                                'Task');
    clear idString RunTime Data Options NeuralSource Sync Task
end

