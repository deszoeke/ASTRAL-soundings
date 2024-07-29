% Parse the university of Wyoming sounding data and put into mat files that
% can be read with Emily's gridding script. These should be named with the
% string <city>_yyyymmddzzzz.mat

% Structures are named "snd" and contain fields with units:
% p (mb), t (C), rh (/100), td (C), ws (knots?), wd (/360), gp (m?)
%% User-defined fields

% Master location of the data
datdir = '/Volumes/cruise/CURRENT/share/data/radiosonde/';

% Names of the cities (folder names) where the radiosonde data is located
citynames = {'Bhubaneswar','Chennai','Kolkata','PortBlair',...
    'Thiruvananthapuram','Vishakhapatnam'};
%% Parse data
% Loop through city name folders. Each file contains one month of data.
% Save the data for each day separately and name with the convention string

% Loop through city folders
for icity = 1:length(citynames)
    disp(['Working on city ' citynames{icity}])
    
    % List the files in the folder
    f = dir([datdir citynames{icity} '/*.txt']);
    f = {f.name};

    % Loop through the files in the folder
    for ifile = 1:length(f)
        disp(['Converting file ' f{ifile}])
        
        % Complete name of the file (with path)
        fname = [datdir citynames{icity} '/' f{ifile}];

        % Read as cell array of strings. each row is a unique string
        dat = table2cell(readtable(fname));
        nrows = length(dat);

        %% Find the indices of the rows that match the header string
        % This tells us where the separation between each balloon is
        headerstr = ['PRES   HGHT   TEMP   DWPT   RELH   MIXR   DRCT'...
            '   SKNT   THTA   THTE   THTV'];
        footerstr = 'Station information and sounding indices';

        iheader = [];
        ifooter = [];
        for irow = 1:nrows
            if strcmpi(dat{irow}, headerstr)
                iheader = cat(2,iheader, irow);
            end

            if strcmpi(dat{irow}, footerstr)
                ifooter = cat(2, ifooter, irow);
            end
        end

        %% Loop through each observation
        % Grab yyyy-mm-dd-zzzz, and the data. Name of the city is already known
        % Indices of the information we want
        % zz   = end-14:end-13
        % day  = end-10:end-9
        % month name = FILENAME(5:6); provided that they're all saved as yyyymm.txt
        % yyyy = end-3:end

        for iday = 1:length(iheader)
            disp(['Saving day ' num2str(iday) '/' num2str(length(iheader))])
            % Grab information about the time so we can give it the correct name to
            % save as
            infostr = dat{iheader(iday) - 2};
            yyyy = infostr(end-3:end);
            mm   = f{ifile}(5:6);
            dd   = infostr(end-10:end-9);
            zzzz = [infostr(end-14:end-13) '00'];

            savename = [datdir 'mat/wyoming/' citynames{icity} '_' yyyy mm dd zzzz];

            % Grab the data. Starts at iheader + 3, and goes until the location of
            % the string "station information and sounding indices"
            subdat = dat(iheader(iday) + 3 : ifooter(iday) - 1);

            % Loop through each string in subdata and grab what we want.
            % Intialize structure fields
            nlevs = length(subdat);
            snd.p  = NaN([1 nlevs]);
            snd.t  = NaN([1 nlevs]);
            snd.rh = NaN([1 nlevs]);
            snd.td = NaN([1 nlevs]);
            snd.ws = NaN([1 nlevs]);
            snd.wd = NaN([1 nlevs]);
            snd.gp = NaN([1 nlevs]);
            for ipres = 1:length(subdat)
                rowsplit = split(subdat{ipres});

                % If the size of the split array is small, this row has a nan value
                % for all fields
                snd.p(ipres) = str2num(rowsplit{1});

                if length(rowsplit) == 11
                    snd.t(ipres)  = str2num(rowsplit{3}); % temperature
                    snd.rh(ipres) = str2num(rowsplit{5}); % relative humidity
                    snd.td(ipres) = str2num(rowsplit{4}); % dewpoint temperature
                    snd.ws(ipres) = str2num(rowsplit{8}); % wind speed
                    snd.wd(ipres) = str2num(rowsplit{7}); % wind direction
                    snd.gp(ipres) = str2num(rowsplit{2}); % height
                end
            end

            % Save snd to the savename and then clear both variables to
            % ensure they're starting clean each time
            save(savename, 'snd')
            clear snd
            clear savename
        end
    end
end