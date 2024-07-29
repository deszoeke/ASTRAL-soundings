% datdir = '/Users/andrea/Documents/miso-bob/soundings/data/';
datadir = '~/Data/cruises/MISOBOB_2019/SR1911/radiosonde/';
casename = 'raw_data_second_20190708_2330.txt';

% Reads formatted text file to a "table" datatype with 16 columns. Each
% column can be accessed by its name. For example:
%   dat.P
% Calling dat.P(1) gives units in cell format. Use curly braces to access
% the actual string: dat.P{1}. The names of each of the columns can be
% accessed with the following: 
%   dat.Properties.VariableNames
dat = readtable([datadir casename]);

VariableNames = dat.Properties.VariableNames;
Units = dat(1:2,:)

% Converting table format to cell format and then to numeric format: an
% example using the pressure column. Using single precision because it
% takes up less space. 
pressure = single(cellfun(@str2num,dat.P(2:end)));