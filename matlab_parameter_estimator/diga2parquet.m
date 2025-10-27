% Converts all .mat files to .parquet (better for loading speed, but files are larger)
relativePath = 'data/validation_new/35'; % path to folder
matFiles = dir(fullfile(relativePath, '*.mat')); % Get all .mat files
outputFolder = 'data/validation_new/35';

for i = 1:length(matFiles)
    filePath = fullfile(relativePath, matFiles(i).name);
    
    % load the .mat file
    load(filePath);
    data = diga.daten; %struct
    field_names = fieldnames(data);
    for j = 1:length(field_names)
        data.(field_names{j}) = data.(field_names{j})';
    end
    
    % save to a relative output folder (set wherever necessary)
    if ~exist(outputFolder, 'dir')
        mkdir(outputFolder); % create folder if it doesn't exist
    end
    
    name = strrep(matFiles(i).name, '.mat', '.parquet');
    parquetwrite(fullfile(outputFolder, name), struct2table(data))

end

clearvars

% OUTPUT: converted .parquet files at 'outputFolder'