% Corrects source for soundings
flist = dir('/Volumes/cruise/SR1911/share/data/radiosonde/mat/*.mat');
for i = 1:length(flist)
    load(fullfile( flist(i).folder, flist(i).name ))
%     snd.source = 'IMD'; % fix in case broken
    fprintf(1,'%s %s\n', flist(i).name, snd.source )
%     save(fullfile( flist(i).folder, flist(i).name ), 'snd')
end

flist = dir('/Volumes/cruise/SR1911/share/data/radiosonde/mat/SallyRide*.mat');
for i = 1:length(flist)
    load(fullfile( flist(i).folder, flist(i).name ))
%     snd.source = 'SR'; % fix in case broken
    fprintf(1,'%s %s\n', flist(i).name, snd.source )
%     save(fullfile( flist(i).folder, flist(i).name ), 'snd')
end

flist = dir('/Volumes/cruise/SR1911/share/data/radiosonde/mat/wyoming/*.mat');
for i = 1:length(flist)
    load(fullfile( flist(i).folder, flist(i).name ))
%     snd.source = 'UWyo'; % fix in case broken
    fprintf(1,'%s %s\n', flist(i).name, snd.source )
%     save(fullfile( flist(i).folder, flist(i).name ), 'snd')
end

