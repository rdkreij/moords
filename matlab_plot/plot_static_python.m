
function plot_static(data, axish, figh, mooring_name, varargin)

% global handles
% programData = getappdata(handles.plot2.fh, 'programData');
% mdcodesOOP = programData.element_database;

% app_level_scale = 1;
% app_level_scale = programData.plot_config.app_level_scale;
% font_name = programData.plot_config.font_name;
% font_size = programData.plot_config.font_size;
% font_weight = programData.plot_config.font_weight;
% app_level_scale = programData.plot_config.app_level_scale;
% rejected_lists = programData.plot_config.rejected_lists;
% rejected_lists = regexp(rejected_lists, ';', 'split');

app_level_scale = 1;
% app_level_scale = 0.5;

cable_scale = 1/10;
cable_scale = 1/15;
cable_scale = 1/20;

font_size = 7;
font_name = 'arial';

include_section_labels = true;
rejected_lists = {};

plot_waves = false;

for ii = 1:2:numel(varargin)
    switch lower(varargin{ii})
        case 'titleaxes'
            titleaxes = varargin{ii+1};
        case 'cable_scale'
            cable_scale = varargin{ii+1};
        case 'rejected_lists'
            rejected_lists = varargin{ii+1};
    end
end

setappdata(figh, 'X', {})
setappdata(figh, 'Y', {})
setappdata(figh, 'Z', {})
setappdata(figh, 'C', {})
setappdata(figh, 'top', {})
setappdata(figh, 'scale', {})
setappdata(figh, 'name', {})
setappdata(figh, 'list', {})

% fh = get_parent_figure(axish);
% setappdata(figh, 'mooring_name', programData.mooring_obj.name)
setappdata(figh, 'mooring_name', mooring_name)
u = uimenu(figh, 'label', 'Save');
u1 = uimenu(u, 'label', 'All formats', 'callback', {@save_mooring_fig, figh});

% include_title = programData.plot_config.include_title;
include_title = true;

if include_title && exist('titleaxes', 'var')
   cla(titleaxes)
   set(titleaxes, 'clipping', 'off', 'visible', 'off')
   args = {'verticalalignment', 'top', 'parent', titleaxes, 'fontweight', 'bold', 'fontname', font_name};
   t = text(0, 1, sprintf('   Experiment: %s', programData.experiment_name), args{:}) ;
   te = get(t, 'extent');
   t = text(0, te(2), sprintf('   Mooring: %s', programData.mooring_obj.name), args{:}) ;
   te = get(t, 'extent');
   t = text(0, te(2), sprintf('   Longitude: %0.4f', programData.mooring_obj.lon), args{:}) ;
   te = get(t, 'extent');
   t = text(0, te(2), sprintf('   Latitude: %0.4f', programData.mooring_obj.lat), args{:}) ;
   te = get(t, 'extent');
   t = text(0, te(2), sprintf('   Depth: %0.0fm', programData.mooring_obj.depth), args{:}) ;
   for ii = 1:numel(programData.mooring_obj.annot.comments_for_header)
       te = get(t, 'extent');
       programData.mooring_obj.annot.comments_for_header{ii} = programData.mooring_obj.annot.comments_for_header{ii};
       comment = programData.mooring_obj.annot.comments_for_header{ii};
       while iscell(comment)
          programData.mooring_obj.annot.comments_for_header{ii} = comment{1};
          comment = comment{1};
       end
       t = text(0, te(2), sprintf('        %s', programData.mooring_obj.annot.comments_for_header{ii} ), args{:}) ;
   end
end
gap = 0.2;
flange = 0.2;
[cylX,cylY,cylZ] = cylinder;
[sphX,sphY,sphZ] = sphere;

start_xl = get(axish, 'xlim');
start_zl = get(axish, 'zlim');

axes(axish)
cla(axish)
hold on

mooring_elements = data.moorele;
clampon_elements = data.mooreleCO;

table = data.bonus_mtab;

[n_rows, n_cols] = size(table);
table = table(n_rows:-1:1, :); % They are backwar

clamp_ind = strcmp(table.("bool_clampon"), 'True');
clamp_table = table(clamp_ind, :);
% clamp_tops = clamp_table.("CLAMP TOP HEIGHT");
clamp_tops = clamp_table.("height");

% mooring_obj = programData.mooring_obj;
% mooring_elements = mooring_obj.mooring_elements;
% clampon_elements = mooring_obj.clampon_elements;

table(1, :);

%% Yuck Solander hacks

% clamp_tops = [cellfun(@str2num, clamp_tops)];

%%
if isempty(mooring_elements)
    return
end

top = 0;
topCO = 0;
mps = [];
% annot = 2;
% COoffset = 0.2;
ec = [0.3 0.3 0.3];
ec = 'none';

%%

if false
    left_label = programData.mooring_obj.annot.left_label;
    left_label_2 = programData.mooring_obj.annot.left_second_label;
    
    right_label = programData.mooring_obj.annot.right_label;
    right_label_2 = programData.mooring_obj.annot.right_second_label;
    
    right_section_label = programData.mooring_obj.annot.right_section_label;
    left_section_label = programData.mooring_obj.annot.left_section_label;
else
    standard_co_offset = 1;

    right_label = 2;
    right_label_2 = 10;
    
    left_label = -right_label - standard_co_offset - 0.5;
    left_label_2 = -right_label_2 - standard_co_offset - 0.5;

    right_section_label = 16;
    left_section_label = -right_section_label - standard_co_offset - 1.5 ;
end
%%

cd = {};
ath = [];
ath2 = [];
athl = [];
athl2 = [];
lph = [];
lphl = [];
toff = 0;
profile on

inline_length_cum = 0;

section_cols = jet; section_cols = section_cols(4:24:end, :) * 0.7;

sections_done = 0;
% total_mooring_height = sum(mooring_obj.inline_lengths(mdcodesOOP)); 
total_mooring_height = sum(data.H(1, :)); 
sech_r = []; sech_l =[];sech_rt = []; sech_lt =[];
step_overs = [];

%% Sort Clampons - needs to be better controlled
% cohASB = programData.mooring_obj.heightsCOasb(programData.element_database);
if false
    cohASB = data.ZCO;
    
    [~, i] = sort(cohASB); %% 
end
% programData.mooring_obj.clampon_elements = programData.mooring_obj.clampon_elements(i);

% cohASB = programData.mooring_obj.heightsCOasb(programData.element_database);
% cohAE = programData.mooring_obj.heightsCOae(programData.element_database);
%%

set(axish, 'units', 'pix')
p_orig = get(axish, 'Position');

inline_length_cum = 0;

% for ii = n_rows:-1:1
for ii = 1:n_rows
    if ii == 10*floor(ii/8);
        drawnow
    end

    entry = table(ii, :);


%     name = mooring_elements(ii).name;
%     list = mooring_elements(ii).list;
%     serial = mooring_elements(ii).serial_number;
%     inv = mooring_elements(ii).inverted;
%     name_for_text = drop_UWA(name);
%     if ~isequal(serial, '?') && ~isempty(serial) && programData.plot_config.include_serial_numbers
%         name_for_text = [name_for_text ' - #' serial];
%     end
%     iscable = mooring_elements(ii).iscable;
%     dbd = pull_db(mooring_elements(ii), mdcodesOOP);
%     drag_form = get_drag_form(mooring_elements(ii), mdcodesOOP);
    
%     scale = app_level_scale * dbd.scale * mooring_elements(ii).scale;

    name = entry.name{1};
    

    list = entry.type{1};
    serial = entry.("serial");
    if iscell(serial)
        serial = serial{1};
    end
    inv = false;
    section_name = entry.("section"){1};

    isclampon = strcmpi(entry.("bool_clampon"){1}, 'True');
    iscable = strcmpi(entry.("bool_line"){1}, 'True');
    
    if isclampon
        continue
    end

    % my_db  = data.bonus_ptab.(list);
    % my_dbi = strcmp(my_db.("name"), name);
    % my_dbe = my_db(my_dbi, :);

    diameter = entry.("diameter");
    width = entry.("width");
    if iscable
        drag_form = 'line';
    else
        if diameter == 0
            drag_form = 'cyl'; % For now
        else
            drag_form = 'sph'; % For now
        end
    end
    scale = 1;

    switch drag_form
        case 'line';
            fc = [0 0 0];
            inline_length = entry.("length");
        case 'cyl';
            fc = [1 0.5 0.1];
            inline_length = entry.("length");
        case 'sph';
            fc = [1 0.5 0.1];
            inline_length = entry.("length");
    end
    
%     inline_length_cum = mooring_obj.cumulative_height(ii-1, mdcodesOOP);

    name = drop_uwa(name);
    name_for_text = name;
    plot_name = name;
    disp(plot_name)

    switch name
        case {'Sercel+frame'}
            plot_name = 'Sercel + frame';
        case {'30in (UWA)' '30in'}
            plot_name = '30in';
        case {'35in (UWA)' '35in'}
            plot_name = '35in';
        case {'150kHz+30in_sph'}
            plot_name = 'rdiqm150+30"dwb';
        case {'75kHz+30in_sph'}
            plot_name = 'rdiqm150+30"dwb';
        case {'75kHz+40in+Ti300kHz'}
            plot_name = 'dwb40"+150khz+300khz';
        case {'Ribuck Dual Kit for Sonardyne ORT'}
            plot_name = 'ribuck dual kit';
        case {'ADV+MP+uSquid+2B'}
            plot_name = 'vector_in_185cm_frame';
        case {'100 cm Signature Frame'}
            plot_name = 'signature_in_100cm_frame';
        case {'3M-WagonW', 'Anchor S245'}
            plot_name = '3ww'; 
        case {'Anchor N280'}
            plot_name = '3ww_ro';
        case {'Anchor W310'}
            plot_name = '2ww_sbs_ro';
        case {'13mmbs-l-16mmds', '13mmBS-L-16mmDS', '16mmBS-L-13mmDS'}
            plot_name = name;
        case {'Lander frame'}
            plot_name = 'lander + sling';
    end

    switch drag_form
        case 'cyl'
            try 
                
                [X, Y, Z, C] = draw_instrument(plot_name, inv);
                fc = ec;
            catch me
                X = 0.5*cylX*inline_length;
                Y = cylY*diameter;
                Z = cylZ*inline_length+top;
                Z = cylZ*inline_length;
                C = Z*0;
                fc = [1, 0, 0];
            end
            dz_plot = inline_length*scale;
            dz_real = inline_length;
            mp = top + dz_real/2;
            tp = top + dz_plot/2;
            
%             ph = plot_static_render(X, Y, Z, C, scale, top, ec, ec, axish, name, list);
%             render_instrument2(X, Y, Z, C, scale, top, ec, fc, gca, name, list);
            render_instrument2(X, Y, Z, C, scale, top, ec, fc, axish, name, list);
           
            po = 2;
            elmtX = diameter*scale;
        case 'sph'
            try 
                [X, Y, Z, C] = draw_instrument(plot_name, inv);
                dz = inline_length*scale;
                dz_plot = inline_length*scale;
                dz_real = inline_length;
            catch me
                X = sphX*(diameter/2)*scale;
                Y = sphY*(diameter/2)*scale;
                Z = sphZ*(inline_length/2)*scale;
                C = Z*0;
                fc = [0.8 0 0];
                dz_plot = inline_length*scale;
                dz_real = inline_length;
            end
            
            mp = top + dz_plot/2;
            tp = mp;
            
%             ph = plot_static_render(X, Y, Z, C, scale, top, ec, fc, axish, name, list);
%           render_instrument2(X, Y, Z, C, scale, top, ec, fc, gca, name, list)
            render_instrument2(X, Y, Z, C, scale, top, ec, fc, axish, name, list);

            po = 2;
            elmtX = 0.5*diameter*scale;
        case 'line'
            
            switch list
                case 'chains'
                    ll = 1;
                    try 
                        [X, Y, Z, C,] = draw_instrument('1m_chain', inv);
                        dz_plot = ll*scale;
                        dz_real = ll;
                        
%                         ph = plot_static_render(X, Y, Z, C, scale, top, ec, ec, axish, name, list);
                        render_instrument2(X, Y, Z, C, scale, top, ec, fc, axish, name, list);

                    catch me
                        hi=1;
                    end
                otherwise
                    ll = min([inline_length 1]); % make it at least 5 m
%                     cable_scale = mooring_elements(ii).cable_scale;
                    ll = max([inline_length*cable_scale ll]); % make it at most l/10 m 
            
                    Z = [0 ll/2-gap/2 nan ll/2 ll/2-gap nan ll/2+gap ll/2 nan ll/2+gap/2 ll]  + top; % plot with spaces
                    X = [0 0 nan -flange flange nan -flange flange nan 0 0];
                    Z = [0 ll]  + top; % plot without spaces
                    X = zeros(size(Z));
                    Y = zeros(size(Z));
            
%                     ph = plot3(X, Y, Z, 'color', fc, 'linewidth', 1.2);
                    W = entry.("width");
%                     W = dbd.dim2;
%                     ph = plot_static_render(cylX*dbd.dim2*scale/100, cylY*dbd.dim2*scale/100, cylZ*ll, [0 0 0], 1, top, [0 0 0], [0 0 0], axish, name, list);
                    render_instrument2(cylX*W*scale/100, cylY*W*scale/100, cylZ*ll, [0 0 0], 1, top, [0 0 0], [0 0 0], axish, name, list);

                    dz_plot = ll;
                    dz_real = ll;
            end
            
            mp = nan;
            tp = top + dz_plot/2;
            
            po = 1;
            elmtX = 0;
        otherwise
            error
    end

    inline_length_cum_prev = inline_length_cum;
    inline_length_cum = inline_length_cum + inline_length;
    
    %% Lable Sections
    section_wrapper = [];
    if include_section_labels && ~isempty(section_name)
        sections_done = sections_done + 1;
        lw = 2; af = 0; 
%         secx = annot*[1 1 1 1]*af-[0.2 0 0 0.2];
        sech_r(end+1) = plot3(right_section_label*[1 1 1 1]-[0.2 0 0 0.2], [0 0 0 0], tp+[-1 -1 1 1]*dz_plot/2, 'color', section_cols(sections_done, :), 'linewidth', lw) ;
        sech_l(end+1) = plot3(left_section_label*[1 1 1 1]+[0.2 0 0 0.2], [0 0 0 0], tp+[-1 -1 1 1]*dz_plot/2, 'color', section_cols(sections_done, :), 'linewidth', lw) ;
        sech_rt(end+1) = text(right_section_label, 0, tp, ['Section ' section_name], 'color', section_cols(sections_done, :), 'rotation', 90, 'Fontsize', font_size, 'fontname', font_name, 'fontweight', 'bold',  'verticalalignment', 'top', 'horizontalalignment', 'center');
        sech_lt(end+1) = text(left_section_label, 0, tp, ['Section ' section_name], 'color', section_cols(sections_done, :), 'rotation', 90, 'Fontsize', font_size, 'fontname', font_name, 'fontweight', 'bold',  'verticalalignment', 'bottom', 'horizontalalignment', 'center');
        section_wrapper = sech_l(end);
        section_wrapper_label = sech_lt(end);
    end
    %%
    
    %% Add inline labels
    if ~any(strcmpi(rejected_lists, list))
        ath(end+1, 1) = text(right_label*app_level_scale, 0, tp,  sprintf('  %s' , name_for_text), 'Fontsize', font_size, 'fontname', font_name, 'verticalalignment', 'middle');
        switch drag_form
            case 'line'
                ath2(end+1) = text(right_label_2*app_level_scale, 0, tp,  sprintf('[%0.1f m sect.]' , inline_length), 'Fontsize', font_size, 'fontname', font_name, 'verticalalignment', 'middle');
            otherwise
%                 ath2(end+1) = text(right_label_2*app_level_scale, 0, tp,  sprintf('[%0.1f m ASB]' , inline_length_cum_prev+dz_real/2), 'Fontsize', font_size, 'fontname', font_name, 'verticalalignment', 'middle');
                ath2(end+1) = text(right_label_2*app_level_scale, 0, tp,  sprintf('[%0.1f to %0.1f m ASB]' , inline_length_cum_prev, inline_length_cum_prev+dz_real), 'Fontsize', font_size, 'fontname', font_name, 'verticalalignment', 'middle');
        end

        
%             spots_r = mooring_obj.annot.get_right_arrow(mooring_obj.mooring_elements(ii).annot.right_arrow_break1, mooring_elements(ii).annot.right_arrow_break2, mooring_elements(ii).annot.right_arrow_break3);
        spots_r = [2, 1.5, 1, 0.5];
        spots_r = spots_r*app_level_scale;
        plot3(spots_r(4), [0], [tp], 'k', 'marker', '<', 'markerfacecolor', 'k', 'markersize', 3)
        lph(end+1, 1) = plot3([spots_r(4) spots_r(3)], [0 0], [tp tp], 'k', 'marker', 'none', 'markerfacecolor', 'k', 'markersize', 2);
        lph(end,   2) = plot3([spots_r(3) spots_r(2)], [0 0], [tp tp], 'k', 'marker', 'none', 'markerfacecolor', 'k', 'markersize', 2);
        lph(end,   3) = plot3([spots_r(2) spots_r(1)], [0 0], [tp tp], 'k', 'marker', 'none', 'markerfacecolor', 'k', 'markersize', 2);
            %%  Add Components
        if false % FIX THIS!!!
            children = has_components(mooring_obj.mooring_elements(ii));
            for mm = 1:numel(children)
               if ~isempty(children(mm).serial_number) && programData.plot_config.include_serial_numbers
                   t_in = sprintf('     %s - #%s', children(mm).name, children(mm).serial_number);
               else
                   t_in = sprintf('     %s' , children(mm).name);
               end
               ath(end+1) = text(right_label*app_level_scale, 0, tp, t_in, 'Fontsize', font_size, 'fontname', font_name, 'color', [0.3 0.3 0.3], 'horizontalalignment', 'left'); 
               ath2(end+1) = text(right_label_2*app_level_scale, 0, tp, ' ');
               lph(end+1, 1) = nan;
               lph(end, 2) = nan;
               lph(end, 3) = nan;
            end
        else
            disp('FIX THIS')
        end
        %%
    end
    %%
    %% Clamp-ons
    

        clamp_line_col = 'k';
                    
        % Find my clamp ons
%         coi = mooring_obj.find_my_clampons(ii, mdcodesOOP);
        coi = clamp_tops > inline_length_cum_prev & clamp_tops <= inline_length_cum;
        ncoi = sum(coi);
        disp([name ' has ' ncoi 'clamp on elements'])

        

    if ncoi >0
        inline_el_start = inline_length_cum_prev;
        sts = 0;

        topCO = topCO + top;
        se_counter = 0;

        heights_this_el_asb = [clamp_tops];
        heights_this_el_ae = heights_this_el_asb - inline_length_cum_prev;
        coi_lin = find(coi)';

        [~, resort_co] = sort(heights_this_el_ae(coi));
        coi_lin = coi_lin(resort_co);

        for kk = coi_lin;
            se_counter = se_counter+1;

            try
                %%
%                 nameCO = clampon_elements(kk).name;
%                 listCO = clampon_elements(kk).list;
%                 name_for_textCO = drop_UWA(nameCO);
%                 iscableCO = clampon_elements(kk).iscable;
%                 dbdCO = pull_db(clampon_elements(kk), mdcodesOOP);
%                 drag_formCO = get_drag_form(clampon_elements(kk), mdcodesOOP);
                entryCO = clamp_table(kk, :);

                nameCO = clamp_table(kk, :).name{1};
                listCO = clamp_table(kk, :).type{1};
                serCO = clamp_table(kk, :).("serial");
                if iscell(serCO)
                    serCO = serCO{1};
                end
                

                iscableCO = strcmpi(entryCO.("bool_line"){1}, 'true');

                % my_dbCO  = data.bonus_ptab.(listCO);
                % my_dbCOi = strcmp(my_dbCO.("name"), nameCO);
                % my_dbCOe = my_dbCO(my_dbCOi, :);

                nameCO = drop_uwa(nameCO);
                name_for_textCO = nameCO;

                diameterCO = clamp_table(kk, :).("diameter");
                widthCO = clamp_table(kk, :).("width");


                if iscableCO
                    drag_formCO = 'line';
                else
                    if diameterCO == 0
                        drag_formCO = 'cyl'; % For now
                    else
                        drag_formCO = 'sph'; % For now
                    end
                end

%                 space_cos_evenly = mooring_obj.mooring_elements(ii).space_cos_evenly;
                space_cos_evenly = true;
                if space_cos_evenly
%                     datum_z = dbdCO.datum_z
                    datum_z = 0;
                    rrl = top - datum_z/100 + dz_plot*se_counter/(sum(coi)+1); % this spaces them appropriatly
                else
    %                 rrl = top + dz_plot*(mooring_obj.clampon_elements(kk).height - inline_el_start)/inline_length; % this spaces them appropriatly;
                    rrl = top + dz_plot*(cohAE(kk) - dbdCO.datum_z/100)/inline_length; % this spaces them appropriatly;
                end
            
%                 scaleCO = dbdCO.scale;
                scaleCO = 1;
%                 scaleCOspec = clampon_elements(kk).scale;
                scaleCOspec = 1;
                scaleCO = app_level_scale * scaleCO * scaleCOspec;
                
                switch drag_formCO
                    case 'line';
                      fc = [1 0.5 0.1];
                      inline_lengthCO = clamp_table(kk, :).("length");
%                     fc = [0 0 0];
%                     inline_lengthCO = clampon_elements(kk).cable_length;
%                     diameterCO = my_dbCOe.("width")/100;
%                     inline_length = my_dbCOe.("width");
%                     error
                case 'cyl';
                    fc = [1 0.5 0.1];
                    inline_lengthCO = clamp_table(kk, :).("length");
                    diameterCO = widthCO;
                case 'sph';
                    fc = [1 0.5 0.1];
                    inline_lengthCO = clamp_table(kk, :).("length");
                end

                plot_name_for_textCO = name_for_textCO;
                switch name_for_textCO
                    case {'SBE39 sync. T'}
                        plot_name_for_textCO = 'sbe39 t';
                end
                
                switch drag_formCO
                    case {'cyl', 'line'}
                        try 
                            [X, Y, Z, C] = draw_instrument(plot_name_for_textCO);
                        catch me
                            X = {0.5*cylX*diameterCO};
                            Y = {cylY*diameterCO};
                            Z = {(cylZ-0.5)*inline_lengthCO};
                            C = {[0.8 0.8 0.8]};
                        end
                        
                        [X, Y, Z] = scale_instrument(X, Y, Z, scaleCO);
%                         [X, Y, Z] = rotate_instrument(X, Y, Z, clampon_elements(kk).dh, clampon_elements(kk).dp, clampon_elements(kk).dr);
                        [X, Y, Z] = shift_instrument(X, Y, Z, -standard_co_offset*scale, 0*scale, 0*scale);

                        render_instrument2(X, Y, Z, C, 1, rrl, ec, ec, axish, name_for_textCO, listCO);
                        po = 2;
                    
                    otherwise
                        keyboard
                end

                
                clear children
                if ~any(strcmpi(rejected_lists, listCO))
%                     spots = mooring_obj.annot.get_left_arrow(mooring_obj.clampon_elements(kk).dx*scale, mooring_obj.clampon_elements(kk).annot.left_arrow_break1, mooring_obj.clampon_elements(kk).annot.left_arrow_break2, mooring_obj.clampon_elements(kk).annot.left_arrow_break3);
%                     spots = spots*app_level_scale;
                    spots = -spots_r - standard_co_offset ;
                    clamp_col = 'k';
%                     if ~programData.plot_config.include_serial_numbers || isempty(mooring_obj.clampon_elements(kk).serial_number)
%                         ns = sprintf(['%s' ], name_for_textCO);
%                     else
%                         ns = sprintf(['#%s - %s' ], mooring_obj.clampon_elements(kk).serial_number, name_for_textCO);
%                     end
                    if lower(nameCO(1))=='s'
                        hi=2;
                    end
                    if isnumeric(serCO)
                        if isnan(serCO)
                            serCO = '';
                        else
                            serCO = num2str(serCO);
                        end
                    end
%                     if isempty(serCO)
                    if ismissing(serCO)
                        ns = sprintf(['%s' ], name_for_textCO);
                    else
                        ns = sprintf(['#%s - %s' ], serCO, name_for_textCO);
                    end
                    athl(end+1) = text(left_label*app_level_scale, 0, rrl, ns, 'Fontsize', font_size, 'fontname', font_name, 'color', clamp_col, 'horizontalalignment', 'right');
% %                     if programData.plot_config.include_height_ae
% %                         t_in = sprintf(['[%0.1f m ASB | %0.1f m AE]' ], cohASB(kk), cohAE(kk));
% %                     else
% %                         t_in = sprintf(['[%0.1f m ASB]' ], cohASB(kk));
% %                     end
                    t_in = sprintf(['[%0.1f m ASB | %0.1f m AE]' ], heights_this_el_asb(kk), heights_this_el_ae(kk));

                    athl2(end+1) = text(left_label_2*app_level_scale, 0, rrl, t_in, 'Fontsize', font_size, 'fontname', font_name, 'color', clamp_col, 'horizontalalignment', 'right');
%                     plot3(spots(4), mooring_obj.clampon_elements(kk).dy, [rrl], 'color', clamp_line_col, 'marker', '>', 'markerfacecolor', clamp_line_col, 'markersize', 3);
                    plot3(spots(4), 0, [rrl], 'color', clamp_line_col, 'marker', '>', 'markerfacecolor', clamp_line_col, 'markersize', 3);

%                     lphl(end+1, 3) = plot3(spots(4:-1:3), mooring_obj.clampon_elements(kk).dy*[1 1], [rrl rrl], 'color', clamp_line_col, 'marker', 'none', 'markerfacecolor', clamp_line_col, 'markersize', 2);
                    lphl(end+1, 3) = plot3(spots(4:-1:3), 0*[1 1], [rrl rrl], 'color', clamp_line_col, 'marker', 'none', 'markerfacecolor', clamp_line_col, 'markersize', 2);

                    lphl(end,   1) = plot3(spots(2:-1:1), [0 0], [rrl rrl], 'color', clamp_line_col, 'marker', 'none', 'markerfacecolor', clamp_line_col, 'markersize', 2);
%                     lphl(end,   2) = plot3(spots(3:-1:2), [mooring_obj.clampon_elements(kk).dy 0], [rrl rrl], 'color', clamp_line_col, 'marker', 'none', 'markerfacecolor', clamp_line_col, 'markersize', 2);
                    lphl(end,   2) = plot3(spots(3:-1:2), [0 0], [rrl rrl], 'color', clamp_line_col, 'marker', 'none', 'markerfacecolor', clamp_line_col, 'markersize', 2);
                
%                     children = has_components(mooring_obj.clampon_elements(kk));
%                     for mm = 1:numel(children)
%                         if ~programData.plot_config.include_serial_numbers || isempty(children(mm).serial_number)
%                             ns = sprintf(['%s' ], children(mm).name);
%                         else
%                             ns = sprintf(['#%s - %s' ], children(mm).serial_number, children(mm).name);
%                         end
%                        athl(end+1) = text(left_label*app_level_scale, 0, rrl, ns, 'Fontsize', font_size, 'fontname', font_name, 'color', [0.3 0.3 0.3], 'horizontalalignment', 'right'); 
%                        athl2(end+1) = text(left_label_2*app_level_scale, 0, rrl, ' ');
%                        lphl(end+1, 1) = nan;
%                        lphl(end, 2) = nan;
%                        lphl(end, 3) = nan;
%                        sts = sts+1;
%                     end
                
                end
                sts = sts+1;
                
            catch me
                disp(me.message)
                for ee = 1:numel(me.stack)
                    disp(me.stack(ee))
                end
                keyboard
            end
        end
        if ~isempty(coi) && ~isempty(section_wrapper) && ~isempty(athl)
           setappdata(section_wrapper, 'sts', sts) 
           setappdata(section_wrapper, 'label', section_wrapper_label) 
           setappdata(athl(end), 'section_wrapper', section_wrapper) 
        end
            
        
        
    %     step_overs(end+1) = step_over;
    else
        disp('FIX THIS')
    end

    top = top + dz_plot;
    mps = cat(1, mps, mp);
%     cd{ii, 1} = ph;
%     cd{ii, 2} = po;
end
profile off
% profile viewer
axis equal off
try % Fails if there are no labels
    ppc_in_y = get(ath(1), 'extent');
catch
end
%% SHIFT LABELS
if isequal(start_xl, [0 1])
    view(0, 0)
end

% zl = [0 top]; set(axish, 'zlim', zl, 'units', 'pix'); p = get(axish, 'Position'); r1 = p(3)/p(4); r = min([r1 0.2]);
% % set(axish, 'zlim', [0 top]*r1/r); 
% xm = 0.5*top/r1; xl = xm*[-1 1]; set(axish, 'xlim', xl); 
% ppc1 = get(ath(1), 'extent');

%%

if plot_waves, % then plot the ocean surface/waves
    depth_ratio = max([total_mooring_height / mooring_obj.depth 0.5]);
    surface_height = top / depth_ratio;
    ampl = 0.003*surface_height;

    xm=get(axish,'XLim');ym=get(axish,'YLim');ym = [-2 2];%ym = xm;
    swx=[min(xm):0.2:max(xm)];
    swy=[min(ym):0.2:max(ym)];
    swx=[min(xm) max(xm)]; dx = min([diff(swx)/10 2]); swx=[min(xm):dx:max(xm)];
    swy=[min(ym) max(ym)]; dy = min([diff(swy)/10 2]); swy=[min(ym):dy:max(ym)];
    
%     box off;
    [wfx, wfy] = meshgrid(swx, swy);
    wfz = ampl*0 * sin(wfx) + ampl*0 * sin(wfy) + ampl*3 * sin(sqrt((5*wfy-200).^2 + (5*wfx-200).^2)/5);
    sh = surf(wfx, wfy, wfz+surface_height, wfz+surface_height, 'parent', axish,'Clipping','off');
    set(sh, 'edgecolor', 'none', 'facecolor', [0 0 1], 'facealpha', 0.1)
    
    sh2 = surf(wfx, wfy, wfz/10, wfz+surface_height, 'parent', axish,'Clipping','off');
    set(sh2, 'edgecolor', 'none', 'facecolor', [0.8 0.6 0], 'facealpha', 1)
    
    zl = [0 surface_height+ampl];
else
    zl = [0 top];
end

%%
r1 = p_orig(3)/p_orig(4); 
xl = [left_section_label-1 right_section_label+1]*app_level_scale;
x_limited = diff(xl)/diff(zl) > r1;
set(axish, 'xlim', xl)
set(axish, 'zlim', zl, 'units', 'pix'); 
% set(axish, 'zlim', [0 top]*r1/r); 
% xm = 0.5*zl(2)*r1; xl = xm*[-1 1]; set(axish, 'xlim', xl); 

dz = 0;



if ~isempty(ath) 
    ppc1 = get(ath(1), 'extent');
    %%
%     ppc = programData.mooring_obj.annot.ppc*ppc1(3)/length(get(ath(1), 'string'));

    yl = get(axish, 'ylim'); set(axish, 'ylim', zl); view(0, 90); ylim(yl)
    
    ppc = get(ath(1), 'extent'); view(0, 0); ppc = ppc(4);
    ppc = 0.5; % Zulberti playing

    for ii = 1:1:numel(ath)-1
        disp([get(ath(ii), 'String') ' vs. ' get(ath(ii+1), 'String')])

        bottom = get(ath(ii), 'position');
        top = get(ath(ii+1), 'position');
        top2 = get(ath2(ii+1), 'position');

        disp(bottom)
    %     top2 = get(ath(ii+3), 'position');
        if bottom(3) + ppc > top(3)
            dz = top(3) - (bottom(3) + ppc);
            set(ath(ii+1), 'position', [top(1:2) top(3) - dz]);
            set(ath2(ii+1), 'position', [top2(1:2) top2(3) - dz]);

            if ~isempty(lph)
                if ~isnan(lph(ii+1, 2))
                    set(lph(ii+1, 2), 'zdata', get(lph(ii+1, 2), 'zdata')+[0 -dz]);
                    set(lph(ii+1, 3), 'zdata', get(lph(ii+1, 3), 'zdata')+[-dz]);
                end
            end

        end
    end
end
dzr = dz;
dz = 0;
for ii = 1:1:numel(athl)
    if ii < numel(athl)
        bottom = get(athl(ii), 'position');
        top = get(athl(ii+1), 'position');
        
        top2 = get(athl2(ii+1), 'position');
        if bottom(3) + ppc > top(3)
            try
                dz = top(3) - (bottom(3) + ppc);
                set(athl(ii+1), 'position', [top(1:2) top(3) - dz]);
                set(athl2(ii+1), 'position', [top2(1:2) top2(3) - dz]);

                set(lphl(ii+1, 2), 'zdata', get(lphl(ii+1, 2), 'zdata')+[0 -dz]);
                set(lphl(ii+1, 1), 'zdata', get(lphl(ii+1, 3), 'zdata')+[-dz]);
            catch
                
            end
        end
    else
%         keyboard
    end
    try
        section_wrapper = getappdata(athl(ii), 'section_wrapper') ;
        zdata = get(section_wrapper, 'zdata'); %zdata(3:4) = zdata(3:4) - dz + ppc/2;
        sts = getappdata(section_wrapper, 'sts');
        section_wrapper_label = getappdata(section_wrapper, 'label'); 
        top = get(athl(ii), 'position');
        bottom = get(athl(ii-sts+1), 'position');
        if bottom(3) > zdata(1)
            zdata(1:2) = bottom(3) - ppc/2;
        end
        if top(3) > zdata(4)
            zdata(3:4) = top(3) + ppc/2;
        end
        set(section_wrapper, 'zdata', zdata)
        p = get(section_wrapper_label, 'position'); p(3) = mean(zdata); set(section_wrapper_label, 'position', p);
    end
end
%%
zlim(get(gca, 'zlim'))
% xlim([-2 1])
i = ~isnan(mps);
[zt, zti] =  unique(mps(~isnan(mps)));

if false
    static_heights = cumsum(mooring_obj.inline_lengths(mdcodesOOP));
    ztl = static_heights(i(end:-1:1));
    ztl = ztl(zti);
    set(axish, 'ztick', zt, 'xtick', -1000:1000:1000, 'ytick', [], 'fontsize', font_size, 'fontname', font_name, 'clipping', 'off')
    set(axish, 'zticklabel', ztl, 'zgrid', 'on')
end


set(axish, 'xlim', xl)
box off
xlabel('Width (mm)', 'fontweight', 'bold', 'fontsize', 9)
zlabel('Height (m ASB)', 'fontweight', 'bold', 'fontsize', 9)

l = light('Position',[0 -250 250]);                 % create a light
lighting gouraud    % prefered method for lighting curved surfaces
material dull 

set(axish, 'userdata', cd)
setappdata(gcf, 'static_texts', ath)
setappdata(gcf, 'static_texts2', ath2)

for tt = 1:numel(ath)
    ext = get(ath(tt), 'extent');
end

zl = get(gca, 'zlim'); 
zl(2) = max([zl(2) 30]); 
set(axish, 'zlim', zl)
% set(axish, 'xlim', diff(zl)*r1*[-0.5 0.5])
set(axish, 'xlim', xl)
% xl = get(gca, 'xlim');
% ratio = max(l)/max(step_overs);
% for tt = 1:numel(sech_rt)
%     p = get(sech_rt(tt), 'Position'); set(sech_rt(tt), 'Position', p+[xl(2) 0 0])
% 	p = get(sech_r(tt), 'xdata'); set(sech_r(tt), 'xdata', p+xl(2))
% end
% for tt = 1:numel(sech_lt)
%     p = get(sech_lt(tt), 'Position'); set(sech_lt(tt), 'Position', p+[xl(1) 0 0])
% 	p = get(sech_l(tt), 'xdata'); set(sech_l(tt), 'xdata', p+xl(1))
% end

if ~isequal(start_xl, [0 1])
    set(axish, 'zlim', start_zl)
    set(axish, 'xlim', start_xl)
end

% X = getappdata(handles.plot2.fh, 'X');
% Y = getappdata(handles.plot2.fh, 'Y');
% Z = getappdata(handles.plot2.fh, 'Z');
% C = getappdata(handles.plot2.fh, 'C');
% scale = getappdata(handles.plot2.fh, 'scale');
% top = getappdata(handles.plot2.fh, 'top');
% name = getappdata(handles.plot2.fh, 'name');
% list = getappdata(handles.plot2.fh, 'list');

if false
    fn = mooring_obj.name;
    save(fullfile(programData.mooring_save_location, latex_legalise2(sprintf('staticImage%s.mat', font_name))), 'X', 'Y', 'Z', 'C', 'top', 'scale', 'name', 'list');
end

% figure;hold on
% for ii = 1:numel(X)
%     plot_static_render(X{ii}, Y{ii}, Z{ii}, C{ii}, scale{ii}, top{ii}, ec, fc, gca)
% end

% error
end

function text = drop_uwa(text)

    try
        text = strip(text);
        if strcmpi(text(length(text)-4:end), '(uwa)')
            text = text(1:length(text)-5);
            text = strip(text);
        end
    end
    try
        text = strip(text);
        if strcmpi(text(length(text)-8:end), '(passive)')
            text = text(1:length(text)-9);
            text = strip(text);
        end
    end
end