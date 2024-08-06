p=pathdef;
path(p);
addpath(['.' filesep]);
addpath(genpath(['.' filesep 'third-party']));
addpath(genpath(['.' filesep 'common']));
addpath(genpath(['.' filesep 'synth_train_generation']));
addpath(genpath(['.' filesep 'data_analysis']));

set(0, 'DefaultAxesFontSize', 20);
set(0, 'DefaultLineLineWidth', 2);
set(0, 'DefaultTextFontSize', 20);
set(0, 'DefaultErrorBarLineWidth', 2);
