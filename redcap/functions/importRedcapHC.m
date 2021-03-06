function dataHC = importRedcapHC(filename, startRow, endRow)
%IMPORTFILE Import numeric data from a text file as a matrix.
%   DATAHC = IMPORTFILE(FILENAME) Reads data from text file FILENAME for
%   the default selection.
%
%   DATAHC = IMPORTFILE(FILENAME, STARTROW, ENDROW) Reads data from rows
%   STARTROW through ENDROW of text file FILENAME.
%
% Example:
%   dataHC = importfile('HCMRINMRStudy_DATA_LABELS_2018-07-09_1335.csv', 2, 8);
%
%    See also TEXTSCAN.

% Auto-generated by MATLAB on 2018/07/09 13:41:03

%% Initialize variables.
delimiter = ',';
if nargin<=2
    startRow = 2;
    endRow = inf;
end

%% Read columns of data as text:
% For more information, see the TEXTSCAN documentation.
formatSpec = '%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%[^\n\r]';

%% Open the text file.
fileID = fopen(filename,'r');

%% Read columns of data according to the format.
% This call is based on the structure of the file used to generate this
% code. If an error occurs for a different file, try regenerating the code
% from the Import Tool.
dataArray = textscan(fileID, formatSpec, endRow(1)-startRow(1)+1, 'Delimiter', delimiter, 'TextType', 'string', 'HeaderLines', startRow(1)-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
for block=2:length(startRow)
    frewind(fileID);
    dataArrayBlock = textscan(fileID, formatSpec, endRow(block)-startRow(block)+1, 'Delimiter', delimiter, 'TextType', 'string', 'HeaderLines', startRow(block)-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
    for col=1:length(dataArray)
        dataArray{col} = [dataArray{col};dataArrayBlock{col}];
    end
end

%% Close the text file.
fclose(fileID);

%% Convert the contents of columns containing numeric text to numbers.
% Replace non-numeric text with NaN.
raw = repmat({''},length(dataArray{1}),length(dataArray)-1);
for col=1:length(dataArray)-1
    raw(1:length(dataArray{col}),col) = mat2cell(dataArray{col}, ones(length(dataArray{col}), 1));
end
numericData = NaN(size(dataArray{1},1),size(dataArray,2));

for col=[1,3,30,34,35,36,37,38,39,40,41,42,43,44,45,46,62,63,64,65,66,67,68,69,70,71]
    % Converts text in the input cell array to numbers. Replaced non-numeric
    % text with NaN.
    rawData = dataArray{col};
    for row=1:size(rawData, 1)
        % Create a regular expression to detect and remove non-numeric prefixes and
        % suffixes.
        regexstr = '(?<prefix>.*?)(?<numbers>([-]*(\d+[\,]*)+[\.]{0,1}\d*[eEdD]{0,1}[-+]*\d*[i]{0,1})|([-]*(\d+[\,]*)*[\.]{1,1}\d+[eEdD]{0,1}[-+]*\d*[i]{0,1}))(?<suffix>.*)';
        try
            result = regexp(rawData(row), regexstr, 'names');
            numbers = result.numbers;
            
            % Detected commas in non-thousand locations.
            invalidThousandsSeparator = false;
            if numbers.contains(',')
                thousandsRegExp = '^\d+?(\,\d{3})*\.{0,1}\d*$';
                if isempty(regexp(numbers, thousandsRegExp, 'once'))
                    numbers = NaN;
                    invalidThousandsSeparator = true;
                end
            end
            % Convert numeric text to numbers.
            if ~invalidThousandsSeparator
                numbers = textscan(char(strrep(numbers, ',', '')), '%f');
                numericData(row, col) = numbers{1};
                raw{row, col} = numbers{1};
            end
        catch
            raw{row, col} = rawData{row};
        end
    end
end

% Convert the contents of columns with dates to MATLAB datetimes using the
% specified date format.
try
    dates{4} = datetime(dataArray{4}, 'Format', 'yyyy-MM-dd', 'InputFormat', 'yyyy-MM-dd');
catch
    try
        % Handle dates surrounded by quotes
        dataArray{4} = cellfun(@(x) x(2:end-1), dataArray{4}, 'UniformOutput', false);
        dates{4} = datetime(dataArray{4}, 'Format', 'yyyy-MM-dd', 'InputFormat', 'yyyy-MM-dd');
    catch
        dates{4} = repmat(datetime([NaN NaN NaN]), size(dataArray{4}));
    end
end

dates = dates(:,4);

%% Split data into numeric and string columns.
rawNumericColumns = raw(:, [1,3,30,34,35,36,37,38,39,40,41,42,43,44,45,46,62,63,64,65,66,67,68,69,70,71]);
rawStringColumns = string(raw(:, [2,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,31,32,33,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,72,73,74,75,76,77]));


%% Replace non-numeric cells with NaN
R = cellfun(@(x) ~isnumeric(x) && ~islogical(x),rawNumericColumns); % Find non-numeric cells
rawNumericColumns(R) = {NaN}; % Replace non-numeric cells

%% Make sure any text containing <undefined> is properly converted to an <undefined> categorical
for catIdx = [1,2,3,4,5,6,7,8,9,10,11,12,14,16,17,19,21,22,23,24,26,27,28,29,30,31,32,40,41,44,50]
    idx = (rawStringColumns(:, catIdx) == "<undefined>");
    rawStringColumns(idx, catIdx) = "";
end

%% Create output variable
dataHC = table;
dataHC.RecordID = cell2mat(rawNumericColumns(:, 1));
dataHC.SubjectID = categorical(rawStringColumns(:, 1));
dataHC.MRIscreeningrecordID = cell2mat(rawNumericColumns(:, 2));
dataHC.Dateofvisit = dates{:, 1};
dataHC.HTN = categorical(rawStringColumns(:, 2));
dataHC.DM = categorical(rawStringColumns(:, 3));
dataHC.CAD = categorical(rawStringColumns(:, 4));
dataHC.Hasthesubjecteverhadapacemakerorimplantablecardiacdevice = categorical(rawStringColumns(:, 5));
dataHC.CHF = categorical(rawStringColumns(:, 6));
dataHC.OnESA = categorical(rawStringColumns(:, 7));
dataHC.OnIViron = categorical(rawStringColumns(:, 8));
dataHC.OnActiveVitaminDorCinacalcet = categorical(rawStringColumns(:, 9));
dataHC.COPD = categorical(rawStringColumns(:, 10));
dataHC.Arrhythmia = categorical(rawStringColumns(:, 11));
dataHC.Connectivetissuedisease = categorical(rawStringColumns(:, 12));
dataHC.Typeofconnectivetissuedisease = rawStringColumns(:, 13);
dataHC.JointDisease = categorical(rawStringColumns(:, 14));
dataHC.Hasthesubjecteverhadajointreplacement = rawStringColumns(:, 15);
dataHC.Hyperlipidemia = categorical(rawStringColumns(:, 16));
dataHC.Liverdisease = categorical(rawStringColumns(:, 17));
dataHC.VarName21 = rawStringColumns(:, 18);
dataHC.Malignancy = categorical(rawStringColumns(:, 19));
dataHC.MalginancyType = rawStringColumns(:, 20);
dataHC.PeripheralVascularDisease = categorical(rawStringColumns(:, 21));
dataHC.Failedrenaltransplant = categorical(rawStringColumns(:, 22));
dataHC.Venousthromboembolicdisease = categorical(rawStringColumns(:, 23));
dataHC.Cerebralvascularevent = categorical(rawStringColumns(:, 24));
dataHC.Other = rawStringColumns(:, 25);
dataHC.Complete = categorical(rawStringColumns(:, 26));
dataHC.Age = cell2mat(rawNumericColumns(:, 3));
dataHC.Gender = categorical(rawStringColumns(:, 27));
dataHC.Race = categorical(rawStringColumns(:, 28));
dataHC.Ethnicity = categorical(rawStringColumns(:, 29));
dataHC.Heightincm = cell2mat(rawNumericColumns(:, 4));
dataHC.InitialWeightkg = cell2mat(rawNumericColumns(:, 5));
dataHC.FinalWeightkg = cell2mat(rawNumericColumns(:, 6));
dataHC.BMI = cell2mat(rawNumericColumns(:, 7));
dataHC.FingerLengthcm = cell2mat(rawNumericColumns(:, 8));
dataHC.FingerCircumferencePREcm = cell2mat(rawNumericColumns(:, 9));
dataHC.FingerCircumferencePOSTcm = cell2mat(rawNumericColumns(:, 10));
dataHC.CalfLengthcm = cell2mat(rawNumericColumns(:, 11));
dataHC.Spacingbetweencalfelectrodescm = cell2mat(rawNumericColumns(:, 12));
dataHC.MajorcalfcircumferencePREcm = cell2mat(rawNumericColumns(:, 13));
dataHC.MajorcalfcircumferencePOSTcm = cell2mat(rawNumericColumns(:, 14));
dataHC.MinorCalfCircumferencePREcm = cell2mat(rawNumericColumns(:, 15));
dataHC.MinorCalfCircumferencePOSTCM = cell2mat(rawNumericColumns(:, 16));
dataHC.WhichsideofthebodywereNMRmeasurementstaken = categorical(rawStringColumns(:, 30));
dataHC.WherewasBPcuffrelativetoNMRmeasurement = categorical(rawStringColumns(:, 31));
dataHC.Whichsideofbodywasbioimpedancetakenon = categorical(rawStringColumns(:, 32));
dataHC.Locationoftherun = rawStringColumns(:, 33);
dataHC.WhereisEasyLog1located = rawStringColumns(:, 34);
dataHC.Whereiseasylog2located = rawStringColumns(:, 35);
dataHC.AnycommentsonOmegatemperatureprobes = rawStringColumns(:, 36);
dataHC.Listeverythingthepersonhadtoeatdrinkexerciselevelsinthe12hoursb = rawStringColumns(:, 37);
dataHC.ListalltripstothebathroomthatthepersonhadbothatCRCandMartinos = rawStringColumns(:, 38);
dataHC.CRCStaff = rawStringColumns(:, 39);
dataHC.Engineer = categorical(rawStringColumns(:, 40));
dataHC.ClinicalStudyStaff = categorical(rawStringColumns(:, 41));
dataHC.Otherpeoplethereduringtherun = rawStringColumns(:, 42);
dataHC.Comments = rawStringColumns(:, 43);
dataHC.Complete1 = categorical(rawStringColumns(:, 44));
dataHC.Sodiumpre = cell2mat(rawNumericColumns(:, 17));
dataHC.BUNpre = cell2mat(rawNumericColumns(:, 18));
dataHC.Creatininepre = cell2mat(rawNumericColumns(:, 19));
dataHC.WBCpre = cell2mat(rawNumericColumns(:, 20));
dataHC.Hemoglobinpre = cell2mat(rawNumericColumns(:, 21));
dataHC.Hematocritpre = cell2mat(rawNumericColumns(:, 22));
dataHC.Plateletspre = cell2mat(rawNumericColumns(:, 23));
dataHC.Osmolalitypre = cell2mat(rawNumericColumns(:, 24));
dataHC.BNPpre = cell2mat(rawNumericColumns(:, 25));
dataHC.Albuminpre = cell2mat(rawNumericColumns(:, 26));
dataHC.Sodiumpost = rawStringColumns(:, 45);
dataHC.BUNpost = rawStringColumns(:, 46);
dataHC.Creatininepost = rawStringColumns(:, 47);
dataHC.Osmolalitypost = rawStringColumns(:, 48);
dataHC.BNPpost = rawStringColumns(:, 49);
dataHC.Complete2 = categorical(rawStringColumns(:, 50));

% For code requiring serial dates (datenum) instead of datetime, uncomment
% the following line(s) below to return the imported dates as datenum(s).

% dataHC.Dateofvisit=datenum(dataHC.Dateofvisit);

