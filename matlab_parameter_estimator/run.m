% PROGRAM START
% Extracts RC params from each SOC range for an entire test profile
% Writes this data in a variable named 'results' and saves it as a .parquet
filePath = 'data/25grad_300';
savePath = 'parametrization/15grad';
files = dir(fullfile(filePath, '*.parquet'));
N_files = length(files);
check = true;

for i = 1:N_files
    disp(['Processing file Nr.', num2str(i)])
    fpath = fullfile(filePath, files(i).name);
    initdata;

    % hold off;
    
    ix_p = [ix_pulses(diff(ix_pulses)>2);...
        ix_pulses(find(diff(ix_pulses)>2)+2)]; % indices of larger drains
    ix_p = sort(ix_p);
    
    % include first and last indices (-20 for safety)
    if ~(isempty(ix_p) && isempty(ix_pulses)) 
        ix_p = [ix_pulses(1); ix_p; ix_pulses(end)] - 20;
    else
        % ix_pulses = (round(numel(vb)/2):numel(vb));
        % ix_p = [round(numel(vb)/2), numel(vb)];
          ix_pulses = (101:round(numel(vb)/2));
        ix_p = [1, round(numel(vb)/2)];
    end

    % reshapes ix_p into a N/2:2 matrix with cosecutive elements in the same row
    ix_p = reshape(ix_p, 2, [])';
    
    
    %~~~~~Testing~~~~~~~~~
    if check
        ix_temp = (ix_pulses(1)-100:ix_pulses(end))';
        plot(t(ix_temp), vb(ix_temp))
        hold on;
        plot(t(ix_temp), soc(ix_temp))
        for j = 1:size(ix_p, 1)
            xregion([t(ix_p(j,1)), t(ix_p(j, 2))]);
        end
        x = input("Constant SOC sections correct? Yes/No [Enter/x]", "s");
        if x == 'x'
            clear x;
            hold off;
            break
        end
        clear x;
        hold off;
        clear ix_temp;
        %~~~~~~~~~~~~~~~~~~~~~
    end
    
    
    %loop over each section
    for j=1:size(ix_p, 1)
        ix=(ix_p(j,1):ix_p(j, 2))';
        plot(t(ix), vb(ix));
        
        optimrc % returns 'xrc' and 'errors'
        calcv;
        
        disp(['xrc=',num2str(xrc)])
        if check
            x = input('OK? Yes/No [Enter/x]', 's');
        else
            x = ' ';
        end
        if x == 'x'
            while true
                disp('Manually select pulses to fit')
                p;
                optimrc;
                b;
                calcv;
                disp(['xrc=',num2str(xrc)])
                y = input('New fit OK? Yes/No [Enter/x]', 's');
                if y=='x'
                    continue    
                else
                    rec;
                    break
                end

            end
            clear y;

        else
            rec %record values, returns 'results'
        end
        clear x;
    end
    clear j;
    clear ix_p;
    
    parquetwrite(fullfile(savePath,files(i).name), results)
    clearvars -except relativePath files N_files savePath check;
    
end

clearvars