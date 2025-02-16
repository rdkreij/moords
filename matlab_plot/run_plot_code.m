%% load mooring configuration and properties dataset
% This script processes and visualizes mooring configurations by loading 
% design data exported from Python, applying specific settings for each mooring, 
% and generating formatted plots with and without rigging. The resulting figures 
% are saved in designated folders. Note that the zlimits (zl) might need
% adjustment per mooring (see also code below).

fprintf('RUNNING PLOT DESIGN\n')
fprintf('Loading design')

% Load .csv files exported from python
% dir = "../output/fieldwork/export/"; 
% mnames = ["W310","N280","S245","L245"]; % file names of exported .csv
% version_name = "fieldwork"; % version name used to save plots

dir = "../output/example/export/"; 
mnames = ["M1","M2"]; % file names of exported .csv
version_name = "example"; % version name used to save plots

msave  = struct(); % allocate structure
for mooring_name = mnames
    msavei = design_to_mdd_python(dir+mooring_name+'.csv');
    msave.(mooring_name) = msavei; % append to structure
end
fprintf(' - succes\n')


fprintf('Start plotting\n')
for mooring_name = mnames

    fprintf('Plotting mooring %s\n',mooring_name)
    switch mooring_name
        case 'W310'
            cable_scale = 1/10;
            mname_l = 'W310 [red]';
        case 'N280'
            cable_scale = 1/10;
            mname_l = 'N280 [blue]';
        case 'S245'
            cable_scale = 1/10;
            mname_l = 'S245 [green]';
        case 'L245'
            cable_scale = 1/5;
            mname_l = 'L245 [yellow]';
        otherwise
            cable_scale = 1/10;
            mname_l = mooring_name;
    end
    if true
        addpath('./plot formatting')
        addpath('./instruments/')
        addpath('./transformations/')
        data = msave.(mooring_name);
    end
    
    % cla(handles.plot2.ah1)
    % programData = getappdata(handles.plot2.fh, 'programData');
    
    % paper_size = programData.plot_config.paper_size;
    % font_name = programData.plot_config.font_name;
    paper_size = 'a4';
    
    switch upper(paper_size)
        case 'A3'
            ph = 42;
            pw = 29.7;
        case 'A4'
            ph = 29.7;
            pw = 21;
    end
    
    % include_title = programData.plot_config.include_title;
    include_title = false;
    if include_title
        fh = figure('color', [1 1 1]);
        dhcm = min([max([0.5 programData.mooring_obj.annot.designated_header_cm]) ph-0.9]);
    %     P = PlotLayout(fh, [0.1 4 0 18 0.5], [0.2 25 0.2]);
        P = PlotLayout(fh, [0 pw 0], [0.2 ph-0.4-dhcm 0 dhcm 0.2]);
        
        ah2 = axes('position', P{1}, 'fontname', font_name);
        ah = axes('position', P{2}, 'fontname', font_name);
        plot_static_python(ah, 'titleaxes', ah2, 'plot_waves', true, 'app_level_scale', 2)

    else
        % Full plot
        fh = figure('color', [1 1 1]);
        P = PlotLayout(fh, [1 pw-2 1], [1 ph-2 1]);
    
        ah = axes('position', P{1});
        plot_static_python(data, ah, fh, 'W310', 'cable_scale', cable_scale)
        ylim(xlim)
        set(gca, 'clipping', 'off')
        switch mooring_name
            case {'W310'}
                zl = [18, 48];
                titlez = 46;
            case {'N280'}
                zl = [18, 55];
                titlez = 40;
            case {'S245'}
                zl = [18, 55];
                titlez = 40;
            otherwise
                zl = [14, 51];
                titlez = 36;
        end
        zlim(zl)
        xl = xlim;
        text(mean(xl), 0, titlez, mname_l, 'fontsize', 10, VerticalAlignment='top', HorizontalAlignment='center')
        fol = ['./' +paper_size +'figures/' +char(version_name) +'/'];
        mkdir(fol)

        if true
    %         print(fh, ['./' paper_size 'figures/' char(mooring_name) '_with_rigging.png'],'-dpng')
            rigging_name = [fol +'/' char(mooring_name) +'_with_rigging.png'];
            print(fh, rigging_name,'-dpng')
        
            % No rigging plot
            fh2 = figure('color', [1 1 1]);
            P = PlotLayout(fh2, [0 pw 0], [0.2 ph-0.4 0.2]);
            ah = axes('position', P{1});
            plot_static_python(data, ah, fh2, 'W310', 'cable_scale', cable_scale, 'rejected_lists', {'rigging'})
            ylim(xlim)
            set(gca, 'clipping', 'off')
            
            zlim(zl)
            xl = xlim;
            text(mean(xl), 0, titlez, mname_l, 'fontsize', 10, VerticalAlignment='top', HorizontalAlignment='center')

            simple_name = [fol +'/' char(mooring_name) +'.png'];
            print(fh2, simple_name,'-dpng')
        end
    end
end

% [~, fn, ~] = fileparts(programData.mooring_obj.file_name);
% if programData.plot_config.save_automatically
%     save_mooring_fig([], [], fh)
% end
