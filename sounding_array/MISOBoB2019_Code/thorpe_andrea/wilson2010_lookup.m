% Function to load the table from Wilson et al. (2010), which gives, for
% the range for n independent identically distributed noise that are 
% normally distributed between 0 and 1, determine the minimum tnr that 
% needs to be surpassed for the overturn. Use the average tnr over the 
% inversion to determine the tnr

function tnr = wilson2010_lookup(p)
    % Returns the tnr that needs to be surpassed for the selected 
    % confidence level
    % Input variables
    %   p - the confidence level (options - 50,90,95,99)
%     datdir = '/Users/andrea/Documents/miso-bob/soundings/papers/WN.txt';
    datdir = '/Users/sdeszoek/Data/cruises/MISOBOB_2019/SR1911/scientific_analysis/programs/radiosonde/thorpe_andrea/WN.txt';
    
    dat = readtable(datdir);
    
    % This is the range of the data. To get the tnr, we also need to divide
    % by n-1, 
    eval(['tnr = dat.W' num2str(p) ';'])
    n = dat.n;
    tnr = tnr./(n-1);

end