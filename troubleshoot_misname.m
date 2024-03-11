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
    hArrayS1 = hst.Array('S1', 'S4');
    hSubjectS4 = hst.Subject('S4', hDebug);
    tempDataFile.NeuralSource.hArrayMaps{1} = hArrayMPF.hArrayMap;
    tempDataFile.NeuralSource.hArrayMaps{2} = hArrayS1.hArrayMap;
    tempDataFile.Options.arrays{1}       = hArrayMPF;
    tempDataFile.Options.arrays{2}       = hArrayS1;
    tempDataFile.Options.subject         = hSubjectS4;
    tempDataFile.Options.output          = replace(tempDataFile.Options.output,'S3','S4'); 
    tempDataFile.NeuralSource.hCBMEX.arrayString = {hArrayMPF.ID, hArrayS1.ID};

    tempDataFile.NeuralSource.output = replace(tempDataFile.NeuralSource.output,'S3','S4');
    tempDataFile.NeuralSource.hCBMEX.outputPath = replace(tempDataFile.NeuralSource.hCBMEX.outputPath,'S3','S4');
    tempDataFile.NeuralSource.hCBMEX.recordFilenames =  replace(tempDataFile.NeuralSource.hCBMEX.recordFilenames,'MPx','MPF');
    tempDataFile.NeuralSource.hCBMEX.recordDirectories = replace(tempDataFile.NeuralSource.hCBMEX.recordDirectories,'MPx','MPF');
    tempDataFile.NeuralSource.hCBMEX.recordFilenames =  replace(tempDataFile.NeuralSource.hCBMEX.recordFilenames,'S3','S4');
    tempDataFile.NeuralSource.hCBMEX.recordDirectories = replace(tempDataFile.NeuralSource.hCBMEX.recordDirectories,'S3','S4');
    for tr = 1:tempDataFile.Task.nTrials
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

