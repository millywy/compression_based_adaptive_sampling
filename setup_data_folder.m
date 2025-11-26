% Setup script to organize data files into a Data/ folder
% This script copies all required .mat files from the TROIKA dataset

clear; clc;

% Define source and destination paths
source_training = './hpi-dhc TROIKA main datasets-IEEE_SPC_2015/Training_data/';
source_test = './hpi-dhc TROIKA main datasets-IEEE_SPC_2015/TestData/';
source_truebpm = './hpi-dhc TROIKA main datasets-IEEE_SPC_2015/TrueBPM/';
dest_folder = './Data/';

% Create Data folder if it doesn't exist
if ~exist(dest_folder, 'dir')
    mkdir(dest_folder);
    fprintf('Created Data folder\n');
end

% Copy training data and BPMtrace files
fprintf('Copying training data files...\n');
training_files = {
    'DATA_01_TYPE01.mat'
    'DATA_02_TYPE02.mat'
    'DATA_03_TYPE02.mat'
    'DATA_04_TYPE02.mat'
    'DATA_05_TYPE02.mat'
    'DATA_06_TYPE02.mat'
    'DATA_07_TYPE02.mat'
    'DATA_08_TYPE02.mat'
    'DATA_09_TYPE02.mat'
    'DATA_10_TYPE02.mat'
    'DATA_11_TYPE02.mat'
    'DATA_12_TYPE02.mat'
};

for i = 1:length(training_files)
    src = [source_training, training_files{i}];
    dst = [dest_folder, training_files{i}];
    if isfile(src)
        copyfile(src, dst);
        fprintf('  Copied %s\n', training_files{i});
    else
        fprintf('  WARNING: %s not found\n', src);
    end
end

% Copy BPMtrace files for training data
fprintf('Copying BPMtrace files...\n');
bpmtrace_files = {
    'DATA_01_TYPE01_BPMtrace.mat'
    'DATA_02_TYPE02_BPMtrace.mat'
    'DATA_03_TYPE02_BPMtrace.mat'
    'DATA_04_TYPE02_BPMtrace.mat'
    'DATA_05_TYPE02_BPMtrace.mat'
    'DATA_06_TYPE02_BPMtrace.mat'
    'DATA_07_TYPE02_BPMtrace.mat'
    'DATA_08_TYPE02_BPMtrace.mat'
    'DATA_09_TYPE02_BPMtrace.mat'
    'DATA_10_TYPE02_BPMtrace.mat'
    'DATA_11_TYPE02_BPMtrace.mat'
    'DATA_12_TYPE02_BPMtrace.mat'
};

for i = 1:length(bpmtrace_files)
    src = [source_training, bpmtrace_files{i}];
    dst = [dest_folder, bpmtrace_files{i}];
    if isfile(src)
        copyfile(src, dst);
        fprintf('  Copied %s\n', bpmtrace_files{i});
    else
        fprintf('  WARNING: %s not found\n', src);
    end
end

% Copy test data files
fprintf('Copying test data files...\n');
test_files = {
    'TEST_S01_T01.mat'
    'TEST_S02_T01.mat'
    'TEST_S02_T02.mat'
    'TEST_S03_T02.mat'
    'TEST_S04_T02.mat'
    'TEST_S05_T02.mat'
    'TEST_S06_T01.mat'
    'TEST_S06_T02.mat'
    'TEST_S07_T02.mat'
    'TEST_S08_T01.mat'
};

for i = 1:length(test_files)
    src = [source_test, test_files{i}];
    dst = [dest_folder, test_files{i}];
    if isfile(src)
        copyfile(src, dst);
        fprintf('  Copied %s\n', test_files{i});
    else
        fprintf('  WARNING: %s not found\n', src);
    end
end

% Copy TrueBPM files (ground truth for test data)
fprintf('Copying TrueBPM files...\n');
truebpm_files = {
    'True_S01_T01.mat'
    'True_S02_T01.mat'
    'True_S02_T02.mat'
    'True_S03_T02.mat'
    'True_S04_T02.mat'
    'True_S05_T02.mat'
    'True_S06_T01.mat'
    'True_S06_T02.mat'
    'True_S07_T02.mat'
    'True_S08_T01.mat'
};

for i = 1:length(truebpm_files)
    src = [source_truebpm, truebpm_files{i}];
    dst = [dest_folder, truebpm_files{i}];
    if isfile(src)
        copyfile(src, dst);
        fprintf('  Copied %s\n', truebpm_files{i});
    else
        fprintf('  WARNING: %s not found\n', src);
    end
end

% Handle Extra_TrainingData
fprintf('Processing Extra_TrainingData.zip...\n');
extra_zip = './hpi-dhc TROIKA main datasets-IEEE_SPC_2015/Extra_TrainingData.zip';
extra_temp = './Extra_TrainingData_temp/';

if isfile(extra_zip)
    % Unzip to temporary folder
    unzip(extra_zip, extra_temp);
    
    % Copy DATA_S04_T01.mat
    src_data = [extra_temp, 'Extra_TrainingData/DATA_S04_T01.mat'];
    dst_data = [dest_folder, 'DATA_S04_T01.mat'];
    if isfile(src_data)
        copyfile(src_data, dst_data);
        fprintf('  Copied DATA_S04_T01.mat\n');
    end
    
    % Copy BPM_S04_T01.mat as DATA_S04_T01_BPMtrace.mat (ground truth)
    src_bpm = [extra_temp, 'Extra_TrainingData/BPM_S04_T01.mat'];
    dst_bpm = [dest_folder, 'DATA_S04_T01_BPMtrace.mat'];
    if isfile(src_bpm)
        copyfile(src_bpm, dst_bpm);
        fprintf('  Copied BPM_S04_T01.mat as DATA_S04_T01_BPMtrace.mat\n');
    end
    
    % Clean up temporary folder
    rmdir(extra_temp, 's');
    fprintf('  Cleaned up temporary files\n');
else
    fprintf('WARNING: Extra_TrainingData.zip not found\n');
end

fprintf('\nData folder setup complete!\n');
fprintf('Total files in Data folder: %d\n', length(dir([dest_folder, '*.mat'])));
