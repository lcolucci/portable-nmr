%% Import data from config.xlsx
% Script for importing data from the following spreadsheet:
%
% Specify the workbook and worksheet. For example: 
%    pathToWorkbook: /Users/linacolucci/Documents/GitHub/MRI/config.xlsx
%    sheetName: MRI_Study_T2
%
% If the structure (# columns, location of columns, etc.) of the excel
% sheet changes, need to update this function accordingly. 
%
% Auto-generated by MATLAB on 2018/03/07 19:32:45 and edited by LAC

function config = load_config_file(pathToWorkbook, sheetName)

%% Checks
if ~ischar(pathToWorkbook) | ~ischar(sheetName)
    error('The filePathToWorkbook or sheetName are not proper character arrays')
end

%% Import the data
[~, ~, raw] = xlsread(pathToWorkbook,sheetName);
raw = raw(2:end,:);
try 
    stringVectors = string(raw(:,[2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20])); %works for R2017a
    stringVectors(ismissing(stringVectors)) = '';
catch 
    stringVectors = raw(:,[2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20]); %works for R2016a
end
raw = raw(:,1);

%% Create output variable
data = reshape([raw{:}],size(raw));

%% Create table
config = table;

%% Allocate imported array to column variable names
config.id = data(:,1);
config.subject = stringVectors(:,1);
config.AMPM = stringVectors(:,2);
config.description = stringVectors(:,3);
config.pathNiftis = stringVectors(:,4);
config.pathMasks = stringVectors(:,5);
config.legShorterTE = stringVectors(:,6);
config.legLongerTE = stringVectors(:,7);
config.fingerShorterTE = stringVectors(:,8);
config.fingerLongerTE = stringVectors(:,9);
config.additionalScan01 = stringVectors(:,10);
config.additionalScan02 = stringVectors(:,11);
config.additionalScan03 = stringVectors(:,12);
config.additionalScan04 = stringVectors(:,13);
config.additionalScan05 = stringVectors(:,14);
config.additionalScan06 = stringVectors(:,15);
config.additionalScan07 = stringVectors(:,16);
config.additionalScan08 = stringVectors(:,17);
config.additionalScan09 = stringVectors(:,18);
config.additionalScan10 = stringVectors(:,19);

%% Clear temporary variables
clearvars data raw stringVectors;