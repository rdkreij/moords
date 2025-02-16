function [ X, Y, Z, C ] = draw_instrument( instrument , inv, varargin)

if nargin < 2
    inv = false;
end
if nargin < 3
    oc = {};
end

[cylX,cylY,cylZ] = cylinder;
[sphX,sphY,sphZ] = sphere;

'                ';

X = {};
Y = {};
Z = {};
C = {};
        
switch lower(instrument)
    
    case {'1m_chain' '20cm_chain'}
        switch lower(instrument)
    
            case '1m_chain' 
                links = 27;
            case '20cm_chain'
                links = 3;
        end
        outer = 0.0485/2;
        inner = 0.008*0.7;
        for ii = 1:links
            if ii/2 ~= round(ii/2)
                [Xo, Yo, Zo] = get_toroid_yz(outer, inner, 12);
            else
                [Xo, Yo, Zo] = get_toroid_xz(outer, inner, 12);
            end
            Zo = Zo + (ii-1)*2*(outer - 1*inner);
            max(Zo(:))
            X = [X Xo];
            Y = [Y Yo];
            Z = [Z Zo];
            C = [C {[0 0 0]}];
        end
         
        return
        
    case {'shac-l-shac'}
        
        outer = 0.04;
        inner = 0.008;
        [Xo, Yo, Zo] = get_toroid_xz(outer, inner, 12);
        X = [X Xo];
        Y = [Y Yo];
        Z = [Z Zo];
        C = [C {[0 0 0]}];
        
        [Xo, Yo, Zo] = get_toroid_yz(outer, inner, 12);
        Zo = Zo + 2*(outer - 1*inner);

        X = [X Xo];
        Y = [Y Yo];
        Z = [Z Zo];
        C = [C {[0 0 0]}];
        
        [Xo, Yo, Zo] = get_toroid_xz(outer, inner, 12);
        Zo = Zo + 4*(outer - 1*inner);
        X = [X Xo];
        Y = [Y Yo];
        Z = [Z Zo];
        C = [C {[0 0 0]}];
        
        return
        [ X, Y, Z, C ] = draw_instrument( 'linkxz' );
        for ii = 1:numel(X)
            Z{ii} = Z{ii} + 0.01;
        end
        [ Xo, Yo, Zo, Co ] = draw_instrument( 'linkyz' );
        for ii = 1:numel(Xo)
            Zo{ii} = Zo{ii} + 0.08 - 2*0.005;
        end
        X = [X Xo];
        Y = [Y Yo];
        Z = [Z Zo];
        C = [C Co];
        [ Xo, Yo, Zo, Co ] = draw_instrument( 'linkxz' );
        for ii = 1:numel(Xo)
            Zo{ii} = Zo{ii} + 2*(0.08 - 2*0.005);
        end
        X = [X Xo];
        Y = [Y Yo];
        Z = [Z Zo];
        C = [C Co];
    
    case {'16mmds-l-16mmds' '16mmbs-l-16mmbs' '13mmbs-l-16mmds' '16mmds-l-13mmbs' '13mmbs-l-13mmbs'}
        
        [x, y, z, c] = draw_instrument('16mmbs', false);
        [X, Y, Z, C] = add_XYZC_to_design(X, Y, Z, C, x, y, z, c);

        [x, y, z, c] = draw_instrument('herc-alloy oblong links 6-10mm', false);
        [x, y, z] = shift_instrument(x, y, z, 0, 0, 0.07);
        [X, Y, Z, C] = add_XYZC_to_design(X, Y, Z, C, x, y, z, c);

        [x, y, z, c] = draw_instrument('16mmbs', true);
        [x, y, z] = shift_instrument(x, y, z, 0, 0, 0.15);
        [X, Y, Z, C] = add_XYZC_to_design(X, Y, Z, C, x, y, z, c);
        
%         [X, Y, Z, C] = add_instrument_to_design(X, Y, Z, C, '16mmbs', 0, 0, 0, 0, 0, 0);
%         [X, Y, Z, C] = add_instrument_to_design(X, Y, Z, C, 'herc-alloy oblong links 6-10mm', 0, 0, 0.07, 0, 0, 0);
%         [X, Y, Z, C] = add_instrument_to_design(X, Y, Z, C, '16mmbs', 0, 0, 0.15, 0, 0, 0, 1);
        
    case {'linkyz'}
        
        angs = 0:20:360;
        r = 0.008;
        o = 0.04;
        for ang = angs;
            rotmat = rotation_matrix(0, ang, 0);
            vector = surf2vector(cylX*r, cylY*r+o, cylZ*3*r);
            vector = rotmat*vector;
            [Xr, Yr, Zr] = vector2surf(vector, 2);
            X{end+1} = Xr;
            Y{end+1} = Yr;
            Z{end+1} = Zr+o;
            C{end+1} = [0 0 0];
        end
        
     case {'linkxz'}
        

        outer = 0.0243;
        inner = 0.0056;
        [Xo, Yo, Zo] = get_toroid_xz(outer, inner, 12);
        X = [X Xo];
        Y = [Y Yo];
        Z = [Z Zo];
        C = [C {[0 0 0]}];
        
    case {'11"' '11in plast-' '30"' '35"' '30in' '35in'}
        
        switch lower(instrument)
            case {'11"' '11in plast-'}
                br = 11*2.54/200;
            case {'30"' '30in'}
                br = 30*2.54/200;
            case {'35"' '35in'}
                br = 35*2.54/200;
        end
            
        X{end+1} = sphX * br;
        Y{end+1} = sphY * br;
        Z{end+1} = sphZ * br + br;
        C{end+1} = [0.8 0.8 0];
        
        X{end+1} = cylX * br;
        Y{end+1} = cylY * br;
        Z{end+1} = cylZ * 0.02 + br - 0.01;
        C{end+1} = [0 0 0];
        
    case {'mclane s0369'}
        
        switch lower(instrument)
            case 'mclane s0369'
                br = 0.76/2;
        end
            
        X{end+1} = sphX * br;
        Y{end+1} = sphY * br;
        Z{end+1} = sphZ * br + br;
        C{end+1} = [0.8 0 0];
        
        [X, Y, Z, C] = add_cyl_to_design_3(X, Y, Z, C, [0 0 br-0.01], [0 0 br-0.01], br+0.02, [0.1 0.1], [1 1] ,[0.8 0 0]);
    
        R1 = 0.06;
        R2 = 0.01;
        [x, y, z] = get_toroid_xz(R1, R2, 20);
        [x, y, z] = shift_instrument({x}, {y}, {z}, 0.4*br, 0, br+0.8*br);
        [X, Y, Z, C] = add_XYZC_to_design(X, Y, Z, C, x, y, z, [0.8 0 0]);
        
        [x, y, z] = get_toroid_xz(R1, R2, 20);
        [x, y, z] = shift_instrument({x}, {y}, {z}, -0.4*br, 0, br-0.8*br-2*R1);
        [X, Y, Z, C] = add_XYZC_to_design(X, Y, Z, C, x, y, z, [0.8 0 0]);
        
        R3 = 0.04;
        [X, Y, Z, C] = add_cyl_to_design_3(X, Y, Z, C, [0 0 -R3], [0 0 2*br+R3], R3, [0.1 0.1], [1 1] ,[0.8 0 0]);

    case {'2ww_sbs'}
        
        n = 2;
        
        h = 0.15;
        lh = 0.05;
        wwc = [0 1 0.9] * 0.7;

        offsets = [-0.4, 0.4];
        
        
        ww = 1;
        tow = ww*h; % Top of wheel

        for xx = 1:2
            off = offsets(xx);

            [X, Y, Z, C] = add_cyl_to_design_3(X, Y, Z, C, [off 0 (ww-1)*h], [off 0 tow], 0.4 , 0, [1 1], wwc);
            [X, Y, Z, C] = add_cyl_to_design_3(X, Y, Z, C, [off 0 tow-lh], [off 0 tow+1e-2], 0.4*1.1 , 0, [1 1], wwc*0.8);

            [X, Y, Z, C] = add_cyl_to_design_3(X, Y, Z, C, [off 0 (ww-1)*h], [off 0 tow+lh], 0.1 , 0, [1 1], wwc);
            %             nn = 9;
%             X{end+1} = sphX(end-nn:end, :)*0.2;
%             Y{end+1} = sphY(end-nn:end, :)*0.2;
%             Z{end+1} = sphZ(end-nn:end, :)*0.05+ww*h;
%             C{end+1} = wwc;
%             
%             X{end+1} = sphX(end-nn:end, :)*0.4;
%             Y{end+1} = sphY(end-nn:end, :)*0.4;
%             Z{end+1} = sphZ(end-nn:end, :)*0.01+ww*h;
%             C{end+1} = wwc;
            
        end

        % Central col
        [X, Y, Z, C] = add_cyl_to_design_3(X, Y, Z, C, [0 0 (1)*h], [0 0 tow+0.1], 0.1 , 0, [1 1], wwc);
        
        % Joiner plate
        min_l = 1.25;
        extend = 0.15;
        [X, Y, Z, C] = add_shs_to_design_3(X, Y, Z, C, [-extend 0 (1)*h], [-extend 0 1*h+0.02], [min_l+2*extend 0.4], 0, [0.95 0.95], wwc*0.5);
%         % Rip-off
%         new_col = [0.8 0.2 0.2];
%         [X, Y, Z, C] = add_shs_to_design_3(X, Y, Z, C, [-0.8 0 (1)*h], [-0.8 0 1*h+0.15], [0.8 0.4], 0, [0.95 0.95], new_col);
%         [ X, Y, Z, C ] = add_instrument_to_design(X, Y, Z, C, '1m_chain' , -0.45, 0, tow+0.15, 0, 0, 0);

        [X, Y, Z, C] = add_cyl_to_design_3(X, Y, Z, C, [0 0 0], [0 0 n*h+0.2], 0.015 , 0, [1 1], [0 0 0]);
        
    case {'3ww_ro'}
        % 2 wagonwheels side by side plus rip off
            
        tow = 3*0.15;
        [ X, Y, Z, C ] = add_instrument_to_design(X, Y, Z, C, '3ww' , 0, 0, 0, 0, 0, 0);
        [ X, Y, Z, C ] = add_instrument_to_design(X, Y, Z, C, 'rip-off' , -0.15, 0, tow, 0, 0, 0);

    case {'2ww_sbs_ro'}
        % 2 wagonwheels side by side plus rip off
            
        tow = 0.15;
        [ X, Y, Z, C ] = add_instrument_to_design(X, Y, Z, C, '2ww_sbs' , 0, 0, 0, 0, 0, 0);
        [ X, Y, Z, C ] = add_instrument_to_design(X, Y, Z, C, 'rip-off' , -0.15, 0, tow, 0, 0, 0);

    case {'rip-off'}

        % Rip-off
        new_col = [0.8 0.2 0.2];
        length = 0.7;
        caseheight = 0.15;
        [X, Y, Z, C] = add_shs_to_design_3(X, Y, Z, C, [-length/2 0 0], [-length/2 0 caseheight], [length 0.4], 0, [0.95 0.95], new_col);
        [ X, Y, Z, C ] = add_instrument_to_design(X, Y, Z, C, '1m_chain' , -0.1, 0, caseheight, 0, 0, 0);

    case {'3ww' '2ww'}
        
        switch lower(instrument)
            case '3ww'
                n = 3;
            case '2ww'
                n = 2;
        end
        h = 0.15;
        lh = 0.05;
        wwc = [0 1 0.9] * 0.7;
%         wwc = [1 0.5 0];
        for ww = 1:n

            [X, Y, Z, C] = add_cyl_to_design_3(X, Y, Z, C, [0 0 (ww-1)*h], [0 0 (ww)*h], 0.4 , 0, [1 1], wwc);
            [X, Y, Z, C] = add_cyl_to_design_3(X, Y, Z, C, [0 0 ww*h-lh], [0 0 ww*h+1e-3], 0.4*1.1 , 0, [1 1], wwc*0.8);
            
            nn = 9;
            X{end+1} = sphX(end-nn:end, :)*0.2;
            Y{end+1} = sphY(end-nn:end, :)*0.2;
            Z{end+1} = sphZ(end-nn:end, :)*0.05+ww*h;
            C{end+1} = wwc;
            
            X{end+1} = sphX(end-nn:end, :)*0.4;
            Y{end+1} = sphY(end-nn:end, :)*0.4;
            Z{end+1} = sphZ(end-nn:end, :)*0.01+ww*h;
            C{end+1} = wwc;
            
        end
        
        [X, Y, Z, C] = add_cyl_to_design_3(X, Y, Z, C, [0 0 0], [0 0 n*h+0.2], 0.015 , 0, [1 1], [0 0 0]);
        
    case '2l1sww'
        
        h = 0.15;
        lh = 0.05;
        wwc = [0 1 0.9] * 0.7;
        wwc = [1 1 0.85] * 0.7;
%         wwc = [1 0.5 0];
        n = 3;
        w = [0.4 0.4 0.25];
        h = 0.15;
        for ww = 1:n

            d = w(ww);
            
            [X, Y, Z, C] = add_cyl_to_design_3(X, Y, Z, C, [0 0 (ww-1)*h], [0 0 (ww)*h], d , 0, [1 1], wwc);
            [X, Y, Z, C] = add_cyl_to_design_3(X, Y, Z, C, [0 0 ww*h-lh], [0 0 ww*h], d*1.1 , 0, [1 1], wwc*0.8);
            
            nn = 9;
            X{end+1} = sphX(end-nn:end, :)*d/2;
            Y{end+1} = sphY(end-nn:end, :)*d/2;
            Z{end+1} = sphZ(end-nn:end, :)*0.05+ww*h;
            C{end+1} = wwc;
            
            X{end+1} = sphX(end-nn:end, :)*d;
            Y{end+1} = sphY(end-nn:end, :)*d;
            Z{end+1} = sphZ(end-nn:end, :)*0.01+ww*h;
            C{end+1} = wwc;
            
        end
        
        d = 0.13;
        [X, Y, Z, C] = add_cyl_to_design_3(X, Y, Z, C, [0 0 0], [0 0 0.5], 0.15 , 0, [1 1], wwc);
        [Xo, Yo, Zo] = get_toroid_xz(0.13, 0.02, 20);
        X{end+1} = Xo;
        Y{end+1} = Yo;
        Z{end+1} = Zo+0.6-2*d;
        C{end+1} = wwc;
    case {'300khzsentinel'}
        
        frame_width = 0.35;
        tab_length = 0.1;
        long_bit = 1 - frame_width - 2*tab_length;
        tilt_ang = 45;
        tilt_length = (frame_width/2)/cosd(tilt_ang);
            
        [ Xo, Yo, Zo, Co ] = draw_instrument( 'sentinel' );
        for ii = 1:numel(Xo)
            Zo{ii} = Zo{ii} + 0.35;
        end
        
        X = [X Xo];
        Y = [Y Yo];
        Z = [Z Zo];
        C = [C Co];
        
        angs = [0 90 180 270];
        for ang = angs;
            
            point0 = [0 0 tab_length];
            point1 = [sind(ang)*frame_width/2 cosd(ang)*frame_width/2 tab_length+frame_width/2];
            point2 = [sind(ang)*frame_width/2 cosd(ang)*frame_width/2 tab_length+long_bit+frame_width/2];
            point3 = [0 0 tab_length+long_bit+frame_width];
            
            [X, Y, Z, C] = add_cyl_to_design_3(X, Y, Z, C, point0, point1, 0.015, 1, [1 1], [0 0 0]);
            [X, Y, Z, C] = add_cyl_to_design_3(X, Y, Z, C, point1, point2, 0.015, 1, [1 1], [0 0 0]);
            [X, Y, Z, C] = add_cyl_to_design_3(X, Y, Z, C, point2, point3, 0.015, 1, [1 1], [0 0 0]);
            
            point2_0 = [0 0 tab_length+long_bit/2+frame_width/2];
            point2_1 = [sind(ang)*frame_width/2 cosd(ang)*frame_width/2 tab_length+long_bit/2+frame_width/2];

            [X, Y, Z, C] = add_cyl_to_design_3(X, Y, Z, C, point2_0, point2_1, 0.015, 1, [1 1], [0 0 0]);
            
        end

        X{end+1} = cylX*0.015;
        Y{end+1} = cylY*0.015;
        Z{end+1} = cylZ*tab_length + 0.9;
        C{end+1} = [0 0 0];
        
        X{end+1} = cylX*0.015;
        Y{end+1} = cylY*0.015;
        Z{end+1} = cylZ*tab_length;
        C{end+1} = [0 0 0];
        
    case {'dualort-sentinel'}
        
        frame_width = 0.35;
        
        [ X, Y, Z, C ] = draw_instrument( 'dualort' );
        [ Xo, Yo, Zo, Co ] = draw_instrument( 'sentinel' );
        for ii = 1:numel(Xo)
            Zo{ii} = Zo{ii} + 0.9;
        end
        X = [X Xo];
        Y = [Y Yo];
        Z = [Z Zo];
        C = [C Co];
        
        angs = [0 180 270];
        for ang = angs;
                       

            rotmat = rotation_matrix(ang, 0, 0);
            vector = surf2vector(cylX*0.015+frame_width/2, cylY*0.015, cylZ*0.78+0.6);
            vector = rotmat*vector;
            [Xr, Yr, Zr] = vector2surf(vector, 2);
            X{end+1} = Xr;
            Y{end+1} = Yr;
            Z{end+1} = Zr;
            C{end+1} = [0 0 0];

            rotmat = rotation_matrix(ang, 0, 135);
            vector = surf2vector(cylX*0.015, cylY*0.015, cylZ * 0.25);
            vector = rotmat*vector;
            [Xr, Yr, Zr] = vector2surf(vector, 2);
            X{end+1} = Xr;
            Y{end+1} = Yr;
            Z{end+1} = Zr + 1.55;
            C{end+1} = [0 0 0];
            
            rotmat = rotation_matrix(ang, 0, 90);
            vector = surf2vector(cylX*0.015, cylY*0.015, cylZ * 0.18);
            vector = rotmat*vector;
            [Xr, Yr, Zr] = vector2surf(vector, 2);
            X{end+1} = Xr;
            Y{end+1} = Yr;
            Z{end+1} = Zr + 1.05;
            C{end+1} = [0 0 0];
        
        end

        X{end+1} = cylX*0.015;
        Y{end+1} = cylY*0.015;
        Z{end+1} = cylZ*0.1 + 1.55;
        C{end+1} = [0 0 0];
        
    case {'ribuck dual kit'}
        
        [X, Y, Z, C] = add_instrument_to_design(X, Y, Z, C, 'dbld 3/8" chain', 0.0, 0, 0, 0, 0, 0);
        [X, Y, Z, C] = add_instrument_to_design(X, Y, Z, C, 'dualort', 0.0, 0, 0.2, 0, 0, 0);
        [X, Y, Z, C] = add_instrument_to_design(X, Y, Z, C, '150mm load ring', 0.0, 0, 0.8, 0, 0, 0);

    case {'dualort'}
        
        [ Xo, Yo, Zo, Co ] = draw_instrument( 'ORT' );
        [Xo, Yo, Zo] = shift_instrument(Xo, Yo, Zo, 0.08, 0, 0);
        [X, Y, Z, C] = add_XYZC_to_design(X, Y, Z, C, Xo, Yo, Zo, Co);
        [Xo, Yo, Zo] = shift_instrument(Xo, Yo, Zo, -2*0.08, 0, 0);
        [X, Y, Z, C] = add_XYZC_to_design(X, Y, Z, C, Xo, Yo, Zo, Co);
        
        [X, Y, Z, C] = add_XYZC_to_design(X, Y, Z, C, Xo, Yo, Zo, Co);
        
        [X, Y, Z, C] = add_shs_to_design_3(X, Y, Z, C, [0 -0.01 0.695], [0 0.01 0.695], [0.05 0.25], 0, [0.95 0.95], [0.8 0.8 0.8]);
        
        [X, Y, Z] = shift_instrument(X, Y, Z, 0, 0, -0.05);
        
    case {'quartermaster', 'longranger'}
        
        falseh = 0.065;
        
        angs = 45+[0:90:270];
        r = 0.1;
        o = sqrt(2)*r+.02;
        su = 0.8;
        for ang = angs;
            n = 9;
            rotmat = rotation_matrix(ang, 0, 20);
            vector = surf2vector(sphX(end-n:end, :)*r+o, sphY(end-n:end, :)*r, sphZ(end-n:end, :)*0.015);
            vector = rotmat*vector;
            [Xr, Yr, Zr] = vector2surf(vector, n+1);
            X{end+1} = Xr;
            Y{end+1} = Yr;
            Z{end+1} = Zr+su - falseh;
            C{end+1} = [0.8 0 0];
            
            vector = surf2vector(cylX*r+o, cylY*r, (cylZ-0.974)*0.07);
            vector = rotmat*vector;
            [Xr, Yr, Zr] = vector2surf(vector, 2);
            X{end+1} = Xr;
            Y{end+1} = Yr;
            Z{end+1} = Zr+su - falseh;
            C{end+1} = [0 0 0.8];
        end
        
        for ii = 1:numel(Z)
            Z{ii} = Z{ii} + 0.3;
        end
        
        [X, Y, Z, C] = add_cyl_to_design_2(X, Y, Z, C, [0 0 0], [0 0 0.9], 0.09 , 0.1 ,[0.95 0.95 0.95]);
        [X, Y, Z, C] = add_cyl_to_design_3(X, Y, Z, C, [0 0 0.9], [0 0 0.9], 0.09 , 1 , [1 1], [0.95 0.95 0.95]);
        [X, Y, Z, C] = add_cyl_to_design_2(X, Y, Z, C, [0 0 0], [0 0 0.02], 0.1 , 0.1 ,[0.2 0.2 1]);
    
    case {'signature 1000 sf'}
        
        [ X, Y, Z, C ] = draw_instrument( 'signature 1000 running' );
        
        falseh = 0.065;
        
        angs = [0:90:270 0];
        r = 0.022;
        os = [0.06*[1 1 1 1] 0];
        su = 0.127;
        
        aziz = [30 30 30 30 0];
        cols = jet;
        cols = cols(4:8:end, :);
        cols=jet(6);
%         colsym=downsample(colsym,4)';
        
        inis = [0.7 0.8 0.9 1 1.1 1.2 ];
        azi = 0;
        o = 0;
        ang = 0;
        stt = 0.25;
        
        for ii = [1 2 3 4 5 6]
            
            ini = inis(ii);
            
            offe = ini*sind(azi);

            startX = (o+offe)*cosd(ang);
            startY = (o+offe)*sind(ang);

            height = ini+0.1;
            offe = height*sind(azi);

            endX = (o+offe)*cosd(ang);
            endY = (o+offe)*sind(ang);

            c = ini*ini*[0.8500 0.3250 0.0980];
            [X, Y, Z, C] = add_cyl_to_design_3(X, Y, Z, C, [startX startY stt+ini], [endX endY height+stt], r*0.8 , 0, [1 1], cols(ii, :));
        
        end
        
    case {'signature 1000 variance'}
        
        [ X, Y, Z, C ] = draw_instrument( 'signature 1000 running' );
        
        falseh = 0.065;
        
        angs = [0:90:270 0];
        r = 0.022;
        os = [0.06*[1 1 1 1] 0];
        su = 0.127;
        
        aziz = [30 30 30 30 0];
        cols = jet;
        cols = cols(4:8:end, :);
        cols=jet(6);
%         colsym=downsample(colsym,4)';
        
        for ii = 1:5;
            
            ang = angs(ii);
            azi = aziz(ii);
            
            o = os(ii);
            
            stt = 0.25;
            
            ini = 1;
            offe = ini*sind(azi);
            
            startX = (o+offe)*cosd(ang);
            startY = (o+offe)*sind(ang);
            
            height = 1.1;
            offe = height*sind(azi);
            
            endX = (o+offe)*cosd(ang);
            endY = (o+offe)*sind(ang);
            
            [X, Y, Z, C] = add_cyl_to_design_3(X, Y, Z, C, [startX startY stt+ini], [endX endY height+stt], r*0.8 , 0, [1 1], cols(4, :));
            
        end
            
    case {'signature 1000 running'}
        
        [ X, Y, Z, C ] = draw_instrument( 'signature 1000' );
        
        falseh = 0.065;
        
        angs = [0:90:270 0];
        r = 0.022;
        os = [0.06*[1 1 1 1] 0];
        su = 0.127;
        
        aziz = [30 30 30 30 0];
        cols = jet;
        cols = cols(4:8:end, :);
        cols=jet(6);
%         colsym=downsample(colsym,4)';
        
        for ii = 1:5;

            ang = angs(ii);
            azi = aziz(ii);
            
            o = os(ii);
            
            stt = 0.25;
            offe = stt*sind(azi);
            blank = 0.1;
            offe = blank*sind(azi);
            
            startX = (o+offe)*cosd(ang);
            startY = (o+offe)*sind(ang);
            
            height = 1.5;
            offe = height*sind(azi);
            
            endX = (o+offe)*cosd(ang);
            endY = (o+offe)*sind(ang);
            
            [X, Y, Z, C] = add_cyl_to_design_3(X, Y, Z, C, [startX startY blank+stt], [endX endY height+stt], r*0.8 , 0, [1 1], [0.8 0.8 0.8]);
            
%             ini = 1;
%             offe = ini*sind(azi);
%             
%             startX = (o+offe)*cosd(ang);
%             startY = (o+offe)*sind(ang);
%             
%             height = 1.1;
%             offe = height*sind(azi);
%             
%             endX = (o+offe)*cosd(ang);
%             endY = (o+offe)*sind(ang);
%             
%             [X, Y, Z, C] = add_cyl_to_design_3(X, Y, Z, C, [startX startY stt+ini], [endX endY height+stt], r*0.8 , 0, [1 1], cols(4, :));
            
        end
        
        
        
    case {'signature 1000'}
        
        [ X, Y, Z, C ] = draw_instrument( 'signaturehat' );
        
        for ii = 1:numel(Z)
            Z{ii} = Z{ii} + 0.2;
        end
        
        X{end+1} = cylX*0.07;
        Y{end+1} = cylY*0.07;
        Z{end+1} = cylZ*0.2;
        C{end+1} = [0.1 0.1 0.1];
    
    case {'signaturehat'}
        
        falseh = 0.065;
        
        X{end+1} = sphX(end-9:end, :)*0.08;
        Y{end+1} = sphY(end-9:end, :)*0.08;
        Z{end+1} = sphZ(end-9:end, :)*0.08*0.5+0.01;
        C{end+1} = [0.1 0.1 0.1];
        
        [X, Y, Z, C] = add_cyl_to_design_3(X, Y, Z, C, [0 0 0], [0 0 0.017], 0.079 , 0, [1 1], [0.1 0.1 0.1]);
        
        angs = 0:90:270;
        r = 0.022;
        o = 0.06;
        su = 0.127;
        for ang = angs;
            n = 9;
            rotmat = rotation_matrix(ang, 0, 26);
            vector = surf2vector(sphX(end-n:end, :)*r+o, sphY(end-n:end, :)*r, sphZ(end-n:end, :)*0.005);
            vector = rotmat*vector;
            [Xr, Yr, Zr] = vector2surf(vector, n+1);
            X{end+1} = Xr;
            Y{end+1} = Yr;
            Z{end+1} = Zr+su - falseh;
            C{end+1} = [1 1 1];
        end
        
        X{end+1} = sphX(end-n:end, :)*r;
        Y{end+1} = sphY(end-n:end, :)*r;
        Z{end+1} = sphZ(end-n:end, :)*0.005+0.047;
        C{end+1} = [1 1 1];
        
    case {'sentinel'}
        
        [ X, Y, Z, C ] = draw_instrument( 'sentinelhat' );
        
        [X, Y, Z] = shift_instrument(X, Y, Z, 0, 0, 0.3);
        
        [X, Y, Z, C] = add_cyl_to_design_3(X, Y, Z, C, [0 0 0], [0 0 0.3], 0.13, 0, [1 1], [0.95 0.95 0.95]);
        [X, Y, Z, C] = add_cyl_to_design_3(X, Y, Z, C, [0 0 0], [0 0 0.02], 0.15, 0, [1 1], [0.2 0.2 1]);
        
    case {'sentinelhat'}
        
        falseh = 0.065;
        
        X{end+1} = sphX(end-3:end, :)*0.3;
        Y{end+1} = sphY(end-3:end, :)*0.3;
        Z{end+1} = sphZ(end-3:end, :)*0.3*1.5-0.325 - falseh;
        C{end+1} = [0.2 0.2 1];
        
        [X, Y, Z, C] = add_cyl_to_design_3(X, Y, Z, C, [0 0 0], [0 0 0.02], 0.15, 0, [1 1], [0.2 0.2 1]);
        
        angs = 0:90:270;
        r = 0.05;
        o = 0.08;
        su = 0.13;
        for ang = angs;
            n = 9;
            rotmat = rotation_matrix(ang, 0, 20);
            vector = surf2vector(sphX(end-n:end, :)*r+o, sphY(end-n:end, :)*r, sphZ(end-n:end, :)*0.01);
            vector = rotmat*vector;
            [Xr, Yr, Zr] = vector2surf(vector, n+1);
            X{end+1} = Xr;
            Y{end+1} = Yr;
            Z{end+1} = Zr+su - falseh;
            C{end+1} = [0.8 0 0];
        end
       
     case {'lisst 200x'}
        
        psl = 0.1;
        
        mc = [0.95 0.95 0.95];
        mc_hivis = [0.5 0.5 0.5];
        
        [X, Y, Z, C] = add_cyl_to_design_2({}, {}, {}, {}, [0 0 0], [0 0 0.65], 0.1 , 0.1 , [0.1 0.1 0.1]);
        
     case {'vector'}
        
        psl = 0.1;
        
        mc = [0.95 0.95 0.95];
        mc_hivis = [0.4 0.4 0.4];
        
        [X, Y, Z, C] = add_cyl_to_design_2({}, {}, {}, {}, [0 0 0.45], [0 0 1.25], 0.05 , 0.1 , mc_hivis);
        [X, Y, Z, C] = add_cyl_to_design_2(X, Y, Z, C, [0 0 0.1], [0 0 0.45], 0.013 , 0.1 , mc_hivis);
        
        [X, Y, Z, C] = add_cyl_to_design_2(X, Y, Z, C, [0 0 0.1], [0 0 0.1]+[psl*cosd(90) psl*sind(90) -0.05], 0.013 , 0.1 , [1 0 0]);
        [X, Y, Z, C] = add_cyl_to_design_2(X, Y, Z, C, [0 0 0.1], [0 0 0.1]+[psl*cosd(210) psl*sind(210) -0.05], 0.013 , 0.1 , [0 0 1]);
        [X, Y, Z, C] = add_cyl_to_design_2(X, Y, Z, C, [0 0 0.1], [0 0 0.1]+[psl*cosd(330) psl*sind(330) -0.05], 0.013 , 0.1 , mc_hivis);
 
%         [X, Y, Z, C] = add_cyl_to_design_2(X, Y, Z, C, [0.25*psl*cosd(90) 0.25*psl*sind(90) 0.1], [0 0 0.1]+[0.75*psl*cosd(90) 0.75*psl*sind(90) -0.75*0.05], 0.013 , 0.1 , [1 0 0]);

    case {'microsquidcluster'}
        
        bcd = 0.025;
        bcl = 0.2;
        
        pl = 0.15;
        pd = 0.008;
        
        sf = 1.5;
        
        [X, Y, Z, C] = add_cyl_to_design_2({}, {}, {}, {}, [-bcd 0 0], [-bcd 0 bcl], bcd , 0.5 , [0.2 0.2 0.2]);
        [X, Y, Z, C] = add_cyl_to_design_2(X, Y, Z, C, [bcd 0 0], [bcd 0 bcl], bcd , 0.5 , [0.2 0.2 0.2]);
        [X, Y, Z, C] = add_cyl_to_design_2(X, Y, Z, C, [0 2*bcd*0.85 0], [0 2*bcd*0.85 bcl], bcd , 0.5 , [0.2 0.2 0.2]);
        
        [X, Y, Z, C] = add_cyl_to_design_2(X, Y, Z, C, [-bcd 0 bcl], [-bcd 0 bcl+pl], pd , 0.5 , [0.8 0.8 0.8]);
%         [X, Y, Z, C] = add_cyl_to_design_2(X, Y, Z, C, [bcd 0 bcl], [-bcd+sf*pd 0 bcl+pl], pd , 0.5 , [0.8 0.8 0.8]);
%         [X, Y, Z, C] = add_cyl_to_design_2(X, Y, Z, C, [0 2*bcd*0.85 bcl], [-bcd+(sf/2)*pd (sf)*pd bcl+pl], pd , 0.5 , [0.8 0.8 0.8]);

    case {'microsquid'}
        
        bcd = 0.025;
        bcl = 0.2;
        
        pl = 0.15;
        pd = 0.008;
        
        sf = 1.5;
        
        [X, Y, Z, C] = add_cyl_to_design_2({}, {}, {}, {}, [-bcd 0 0], [-bcd 0 bcl], bcd , 0.5 , [0.2 0.2 0.2]);
%         [X, Y, Z, C] = add_cyl_to_design_2(X, Y, Z, C, [bcd 0 0], [bcd 0 bcl], bcd , 0.5 , [0.2 0.2 0.2]);
%         [X, Y, Z, C] = add_cyl_to_design_2(X, Y, Z, C, [0 2*bcd*0.85 0], [0 2*bcd*0.85 bcl], bcd , 0.5 , [0.2 0.2 0.2]);
        
        [X, Y, Z, C] = add_cyl_to_design_2(X, Y, Z, C, [-bcd 0 bcl], [-bcd 0 bcl+pl], pd , 0.5 , [0.8 0.8 0.8]);
        
    case {'signature_in_100cm_frame'}

        [X, Y, Z, C] = add_instrument_to_design(X, Y, Z, C, '100cm adv frame only', 0.02, 0, 0.06, 0, 0, 0);
        [X, Y, Z, C] = add_instrument_to_design(X, Y, Z, C, 'signature 1000', 0, 0, 0.7, 0, 0, 0);

    case {'vector_in_185cm_frame'}
        
        [X, Y, Z, C] = add_instrument_to_design(X, Y, Z, C, '185cm adv frame only', 0.02, 0, 0.06, 0, 0, 0);
        [X, Y, Z, C] = add_instrument_to_design(X, Y, Z, C, 'microsquid', 0.02, 0, 0.16, 0, 0, 0);
        [X, Y, Z, C] = add_instrument_to_design(X, Y, Z, C, 'vector', 0, 0, 0.5, 0, 0, 0);
        
    case {'vector_in_165cm_frame'}
        
        [X, Y, Z, C] = add_instrument_to_design(X, Y, Z, C, '165cm adv frame only', 0.02, 0, 0.06, 0, 0, 0);
        [X, Y, Z, C] = add_instrument_to_design(X, Y, Z, C, 'microsquid', 0.02, 0, 0.06, 0, 0, 0);
        [X, Y, Z, C] = add_instrument_to_design(X, Y, Z, C, 'vector', 0, 0, 0.4, 0, 0, 0);
    
    case {'185cm adv frame only' '165cm adv frame only' '100cm adv frame only'}
        
        switch lower(instrument)
            case '185cm adv frame only'
                fl = 1.85;
                mp1 = fl - 0.75;
                mp2 = fl - 0.5;
            case '165cm adv frame only'
                fl = 1.65;
                mp1 = fl - 0.75;
                mp2 = fl - 0.5;
            case '100cm adv frame only'
                fl = 1.00;
                mp1 = fl - 0.3;
                mp2 = fl - 0.55;
        end
        
        fw = 0.50;
        ect = 0.05;
        
        fpt = 0.04;
        fpew = 0.0;
        pt = 0.015;
        ro = 0.04;
        
        
        
        mc  = [0.4 0.4 0.4];
        
        [X, Y, Z, C] = add_cyl_to_design_3(X, Y, Z, C, [-fw/2+pt/2+ro -fw/2+pt/2+ro 0], [-fw/2+pt/2+ro -fw/2+pt/2+ro fl], pt, 0, [1 1], mc);
        [X, Y, Z, C] = add_cyl_to_design_3(X, Y, Z, C, [fw/2-pt/2-ro -fw/2+pt/2+ro 0], [fw/2-pt/2-ro -fw/2+pt/2+ro fl], pt, 0, [1 1], mc);
        [X, Y, Z, C] = add_cyl_to_design_3(X, Y, Z, C, [-fw/2+pt/2+ro -(-fw/2+pt/2+ro) 0], [-fw/2+pt/2+ro -(-fw/2+pt/2+ro) fl], pt, 0, [1 1], mc);
        [X, Y, Z, C] = add_cyl_to_design_3(X, Y, Z, C, [fw/2-pt/2-ro -(-fw/2+pt/2+ro) 0], [fw/2-pt/2-ro -(-fw/2+pt/2+ro) fl], pt, 0, [1 1], mc);
        
        [X, Y, Z, C] = add_shs_to_design_3(X, Y, Z, C, [0 0 0], [0 0 ect], fw*[1 1], [0.1 0], [1 1], mc);
        [X, Y, Z, C] = add_shs_to_design_3(X, Y, Z, C, [0 0 fl-ect], [0 0 fl], fw*[1 1], [0 0.1], [1 1], mc);
        
        [X, Y, Z, C] = add_shs_to_design_3(X, Y, Z, C, [0 0 mp1], [0 0 mp1+fpt], fw+fpew*[1 1], 0, [1 1], [0.2 0.2 0.2]);
        [X, Y, Z, C] = add_shs_to_design_3(X, Y, Z, C, [0 0 mp2], [0 0 mp2+fpt], fw+fpew*[1 1], 0, [1 1], [0.2 0.2 0.2]);
        
    case 'ort'
        
        ortl = 0.695;
        caph = 0.14/2;
        capw = 0.14/2;
        long_bit = 0.695 - 2*caph;
        crossw = 0.005;
        crossl = caph/sind(45);%-crossw;
        crossl = caph;%-crossw;
        tab_length = 0;
        col = [0.8 0.8 0.8];col = [0.0 0.0 0.0];
        
        main_tube = (long_bit - 1.5*caph);
        
        angs = 0:360/6:360;
        for ang = angs;
            
            point0 = [0 0 tab_length];
            point1 = [sind(ang)*crossl cosd(ang)*crossl tab_length+caph];
            point2 = [sind(ang)*crossl cosd(ang)*crossl tab_length+long_bit+caph];
            point3 = [0 0 tab_length+long_bit+caph*2];
            
            [X, Y, Z, C] = add_cyl_to_design_3(X, Y, Z, C, point1, point2, crossw, 1, [1 1], col);
            [X, Y, Z, C] = add_cyl_to_design_3(X, Y, Z, C, point2, point3, crossw, 1, [1 1], col);
           
            
        end
       
        [X, Y, Z, C] = add_cyl_to_design_3(X, Y, Z, C, [0 0 caph], [0 0 caph+2*crossw], capw *1.12, 0, [1 1], [197 179 88]/255);
        [X, Y, Z, C] = add_cyl_to_design_3(X, Y, Z, C, [0 0 main_tube+caph], [0 0 main_tube+caph+2*crossw], capw *1.12, 0, [1 1], [197 179 88]/255);
        
        [X, Y, Z, C] = add_cyl_to_design_3(X, Y, Z, C, [0 0 main_tube+caph], [0 0 main_tube+caph+caph], capw * 0.5, 0.6, [1 1], [0 0 0]);
        [X, Y, Z, C] = add_cyl_to_design_3(X, Y, Z, C, [0 0 0], [0 0 caph], capw * 0.6, 0, [1 1], [0.4 0.4 0.4]);

        [X, Y, Z, C] = add_cyl_to_design_3(X, Y, Z, C, [0 0 caph], [0 0 main_tube+caph], (capw - crossw), 0, [1 1], [0.99 0.99 0.99]);

    case {'dwb40"+150khz+300khz'}

    
        [X, Y, Z, C] = add_instrument_to_design(X, Y, Z, C, 'dwb40"+150khz', 0, 0, 0, 0, 0, 0);
        [X, Y, Z, C] = add_inverted_instrument_to_design(X, Y, Z, C, 'sentinel', 0, 0, 0.3, 0, 0, 0);

    case {'rdiqm150+36"dwb' 'rdiqm150+30"dwb' 'dwb40"+150khz'}
        
        switch lower(instrument)
            case 'dwb40"+150khz'
                buoycol = [1 0.5 0];
                br = 1.0160/2;
                frame_length = 1.80;
                adcpshift = 0.5;
                [ Xo, Yo, Zo, Co ] = draw_instrument( 'quartermaster' );
                frame = '1.8m adcp buoy frame';
            case 'rdiqm150+36"dwb'
                buoycol = [1 0.5 0];
                br = .914400/2;
                frame_length = 1.63;
                adcpshift = 0.35;
                [ Xo, Yo, Zo, Co ] = draw_instrument( 'quartermaster' );
                frame = '1.6m adcp buoy frame';
            case 'rdiqm150+30"dwb'
                buoycol = [0.8 0.8 0];
                br = .76/2;
                frame_length = 1.63;
                adcpshift = 0.35;
                [ Xo, Yo, Zo, Co ] = draw_instrument( 'quartermaster' );
                frame = '1.6m adcp buoy frame';
            otherwise
                buoycol = [0.8 0.8 0];
                br = .76/2;
                frame_length = 1.63;
                adcpshift = 0.97;
                [ Xo, Yo, Zo, Co ] = draw_instrument( 'sentinel' );
                frame = '1.6m adcp buoy frame';
        end
        [ X, Y, Z, C ] = draw_instrument( frame );
        
        for ii = 1:numel(Xo)
            Zo{ii} = Zo{ii} + adcpshift;
        end
        X = [X Xo];
        Y = [Y Yo];
        Z = [Z Zo];
        C = [C Co];
        
        X{end+1} = sphX(3:end-2, :) * br;
        Y{end+1} = sphY(3:end-2, :) * br;
        Z{end+1} = sphZ(3:end-2, :) * br + frame_length*0.5;
        C{end+1} = buoycol;
        
        X{end+1} = cylX * br;
        Y{end+1} = cylY * br;
        Z{end+1} = cylZ * 0.02 + frame_length*0.5 - 0.01;
        C{end+1} = [0 0 0.8];
        
    case {'1.6m adcp buoy frame' '1.8m adcp buoy frame'}
        
        switch lower(instrument)
            case '1.6m adcp buoy frame'
                frame_length = 1.64;
            case '1.8m adcp buoy frame'
                frame_length = 1.8;
        end
        caph = 0.45/2;
        long_bit = frame_length - 2*caph;
        crossw = 0.01;
        crossl = caph/sind(45);%-crossw;
        tab_length = 0;
        col = [0.5 0.5 0.5];
        
        angs = [0 90 180 270];
        for ang = angs;
            
            point0 = [0 0 tab_length];
            point1 = [sind(ang)*caph cosd(ang)*caph tab_length+caph];
            point2 = [sind(ang)*caph cosd(ang)*caph tab_length+long_bit+caph];
            point3 = [0 0 tab_length+long_bit+caph*2];
            
            [X, Y, Z, C] = add_cyl_to_design_3(X, Y, Z, C, point0, point1, crossw, 1, [1 1], col);
            [X, Y, Z, C] = add_cyl_to_design_3(X, Y, Z, C, point1.*[0 0 1], point1, crossw, 1, [1 1], col);
            [X, Y, Z, C] = add_cyl_to_design_3(X, Y, Z, C, point1, point2, crossw, 1, [1 1], col);
            [X, Y, Z, C] = add_cyl_to_design_3(X, Y, Z, C, point2.*[0 0 1], point2, crossw, 1, [1 1], col);
            [X, Y, Z, C] = add_cyl_to_design_3(X, Y, Z, C, point2, point3, crossw, 1, [1 1], col);
            
            point2_0 = [0 0 tab_length+long_bit/2+crossl];
            point2_1 = [sind(ang)*crossl cosd(ang)*crossl tab_length+long_bit/2+crossl];

            [X, Y, Z, C] = add_cyl_to_design_3(X, Y, Z, C, point2_0, point2_1, crossw, 1, [1 1], col);
            
        end
    
    case {'sercel + frame'}
            
        [X, Y, Z, C] = add_instrument_to_design(X, Y, Z, C, 'sercel frame', 0, 0, 0, 0, 0, 0);
        [X, Y, Z, C] = add_instrument_to_design(X, Y, Z, C, 'sercel beacon', 0, 0, 0.15, 0, 0, 0);
        
    case {'sercel frame'}
        
        frame_length = 0.5;
        frame_width = 0.1;
        tab_length = 0.03;
        mc  = [0.4 0.4 0.4];
        
        pw_x = 0.02;
        pw_y = 0.01;
        
        long_bit = frame_length - 2*tab_length;
        
        [X, Y, Z, C] = add_shs_to_design_3(X, Y, Z, C, [0 0 0], [0 0 tab_length], [pw_x pw_y] , 0, [0.95 0.95], mc);
        [X, Y, Z, C] = add_shs_to_design_3(X, Y, Z, C, [-frame_width/2 0 tab_length], [-frame_width/2 0 tab_length+long_bit], [pw_x pw_y] , 0, [0.95 0.95], mc);
        [X, Y, Z, C] = add_shs_to_design_3(X, Y, Z, C, [frame_width/2 0 tab_length], [frame_width/2 0 tab_length+long_bit], [pw_x pw_y] , 0, [0.95 0.95], mc);
        [X, Y, Z, C] = add_shs_to_design_3(X, Y, Z, C, [0 0 tab_length+long_bit], [0 0 frame_length], [pw_x pw_y] , 0, [0.95 0.95], mc);
        
        for percs = 0:25:100
            [X, Y, Z, C] = add_shs_to_design_3(X, Y, Z, C, [-(frame_width+pw_x)/2 0 tab_length+percs*long_bit/100], [(frame_width+pw_x)/2 0 tab_length+percs*long_bit/100], [pw_x pw_y] , 0, [0.95 0.95], mc);
        end
        
    case 'lander'
        
        bw = 2.4;
        h = 1.3;
        sink = -0.2;
        d = 0.035;
        lc = [0 0.2 0.6];
        
        rr = 2;
        xshi = (bw - bw*(1/rr))/(1+1/cosd(60));
        zshi = xshi/cosd(60);
        
        sbw = 2.4/rr;
        blo = (bw - sbw)/2;
        
        lt = [0 bw bw/2 ; 0 0 -bw*cosd(30) ; sink sink sink];
        st = [blo blo+bw/2 bw/2 ; -xshi -xshi -xshi-sbw*cosd(30); sink sink sink];
        sth = st + [0 0 0; 0 0 0; h h h];
        
        %%
        [X, Y, Z, C] = add_cyl_to_design(X, Y, Z, C, lt(:, 1), lt(:, 2), d ,lc);
        [X, Y, Z, C] = add_cyl_to_design(X, Y, Z, C, lt(:, 2), lt(:, 3), d ,lc);
        [X, Y, Z, C] = add_cyl_to_design(X, Y, Z, C, lt(:, 3), lt(:, 1), d ,lc);
        %%
        [X, Y, Z, C] = add_cyl_to_design(X, Y, Z, C, st(:, 1), st(:, 2), d ,lc);
        [X, Y, Z, C] = add_cyl_to_design(X, Y, Z, C, st(:, 2), st(:, 3), d ,lc);
        [X, Y, Z, C] = add_cyl_to_design(X, Y, Z, C, st(:, 3), st(:, 1), d ,lc);
        
        %%
        [X, Y, Z, C] = add_cyl_to_design(X, Y, Z, C, sth(:, 1), sth(:, 2), d ,lc);
        [X, Y, Z, C] = add_cyl_to_design(X, Y, Z, C, sth(:, 2), sth(:, 3), d ,lc);
        [X, Y, Z, C] = add_cyl_to_design(X, Y, Z, C, sth(:, 3), sth(:, 1), d ,lc);

        %%
        [X, Y, Z, C] = add_cyl_to_design(X, Y, Z, C, st(:, 1), st(:, 1)+[0 0 h]', d ,lc);
        [X, Y, Z, C] = add_cyl_to_design(X, Y, Z, C, st(:, 2), st(:, 2)+[0 0 h]', d ,lc);
        [X, Y, Z, C] = add_cyl_to_design(X, Y, Z, C, st(:, 3), st(:, 3)+[0 0 h]', d ,lc);
        
        %%
        [X, Y, Z, C] = add_cyl_to_design(X, Y, Z, C, lt(:, 1), st(:, 1)+[0 0 h]', d ,lc);
        [X, Y, Z, C] = add_cyl_to_design(X, Y, Z, C, lt(:, 2), st(:, 2)+[0 0 h]', d ,lc);
        [X, Y, Z, C] = add_cyl_to_design(X, Y, Z, C, lt(:, 3), st(:, 3)+[0 0 h]', d ,lc);
        %% BRACES
        fa = 0.25;
        p1 = fa*(lt(:, 2) - lt(:, 1))+lt(:, 1);
        [X, Y, Z, C] = add_cyl_to_design(X, Y, Z, C, p1, st(:, 1), d ,lc);
        p1 = (1 - fa)*(lt(:, 2) - lt(:, 1))+lt(:, 1);
        [X, Y, Z, C] = add_cyl_to_design(X, Y, Z, C, p1, st(:, 2), d ,lc);
        p1 = fa*(lt(:, 3) - lt(:, 2))+lt(:, 2);
        [X, Y, Z, C] = add_cyl_to_design(X, Y, Z, C, p1, st(:, 2), d ,lc);
        p1 = (1 - fa)*(lt(:, 3) - lt(:, 2))+lt(:, 2);
        [X, Y, Z, C] = add_cyl_to_design(X, Y, Z, C, p1, st(:, 3), d ,lc);
        p1 = fa*(lt(:, 1) - lt(:, 3))+lt(:, 3);
        [X, Y, Z, C] = add_cyl_to_design(X, Y, Z, C, p1, st(:, 3), d ,lc);
        p1 = (1 - fa)*(lt(:, 1) - lt(:, 3))+lt(:, 3);
        [X, Y, Z, C] = add_cyl_to_design(X, Y, Z, C, p1, st(:, 1), d ,lc);
        
        [X, Y, Z] = shift_instrument(X, Y, Z, -mean(lt(1, :)), -mean(lt(2, :)), 0);
        [X, Y, Z] = rotate_instrument(X, Y, Z, 180, 0, 0);
        
    case 'lander sling'
        
        bw = 2.4;
        h = 0;
        sink = -0.2;
        
        rr = 2;
        xshi = (bw - bw*(1/rr))/(1+1/cosd(60));
        
        sbw = bw/rr;
        blo = (bw - sbw)/2;
        
        lt = [0 bw bw/2 ; 0 0 -bw*cosd(30) ; sink sink sink];
        st = [blo blo+bw/2 bw/2 ; -xshi -xshi -xshi-sbw*cosd(30); sink sink sink];
        sth = st + [0 0 0; 0 0 0; h h h];
        
        [X, Y, Z, C] = add_cyl_to_design_2(X, Y, Z, C, sth(:, 1), [mean(sth(1, :)) mean(sth(2, :)) h+1.4]', 0.01 , 0, [0 0 0]);
        [X, Y, Z, C] = add_cyl_to_design_2(X, Y, Z, C, sth(:, 2), [mean(sth(1, :)) mean(sth(2, :)) h+1.4]', 0.01 , 0, [0 0 0]);
        [X, Y, Z, C] = add_cyl_to_design_2(X, Y, Z, C, sth(:, 3), [mean(sth(1, :)) mean(sth(2, :)) h+1.4]', 0.01 , 0, [0 0 0]);
        
        [X, Y, Z] = shift_instrument(X, Y, Z, -mean(lt(1, :)), -mean(lt(2, :)), 0);
        [X, Y, Z] = rotate_instrument(X, Y, Z, 180, 0, 0);
        
    case 'prelude lander'
        
        [ X, Y, Z, C ] = draw_instrument( 'lander' );
%         [X, Y, Z, C] = add_instrument_to_design(X, Y, Z, C, 'quartermaster', -1, 0.4, 0, 0, 0, 0);
%         [X, Y, Z, C] = add_instrument_to_design(X, Y, Z, C, 'vector', -0.5, 0.35, 1.05, 0, 0, 90, true);
        
%         [X, Y, Z] = shift_instrument(X, Y, Z, -1, 0.5, 0);

    case 'lander + sling'
            
        [ X, Y, Z, C ] = draw_instrument( 'lander' );
        [X, Y, Z, C] = add_instrument_to_design(X, Y, Z, C, 'lander sling', 0, 0, 1.1, 0, 0, 0);
    %         [X, Y, Z, C] = add_instrument_to_design(X, Y, Z, C, 'quartermaster', -1, 0.4, 0, 0, 0, 0);
    %         [X, Y, Z, C] = add_instrument_to_design(X, Y, Z, C, 'vector', -0.5, 0.35, 1.05, 0, 0, 90, true);
            
    %         [X, Y, Z] = shift_instrument(X, Y, Z, -1, 0.5, 0);


    case 'prelude lander sling'
        
        [ X, Y, Z, C ] = draw_instrument( 'lander sling' );
%         [X, Y, Z, C] = add_instrument_to_design(X, Y, Z, C, 'quartermaster', -1, 0.4, 0, 0, 0, 0);
%         [X, Y, Z, C] = add_instrument_to_design(X, Y, Z, C, 'vector', -0.5, 0.35, 1.05, 0, 0, 90, true);
        
%         [X, Y, Z] = shift_instrument(X, Y, Z, -1, 0.5, 0);
    
    case {'edgetech sport pop-up' 'pop-up buoy' 'popup buoy'}
        
        [X, Y, Z, C] = add_cyl_to_design_3(X, Y, Z, C, [0 0 0], [0 0 0.7], 0.033 , 0, [0.85 0.95], [0 0 0]);
        [X, Y, Z, C] = add_cyl_to_design_3(X, Y, Z, C, [0 0 0.1], [0 0 0.375], 0.19 , 0, [0.85 0.95], [1 0.8 0]);
        [X, Y, Z, C] = add_cyl_to_design_3(X, Y, Z, C, [0 0 0.375], [0 0 0.66], 0.19 , 0, [1 0.85], [01 0.8 0]);

    case {'16mmbs' '16mmds'}
        
        switch lower(instrument)
            case {'16mmbs' '16mmds'}
                outer = 0.035;
                inner = 0.01;
        end
        
        [Xo, Yo, Zo] = get_toroid_xz(outer, inner, 13);
        Xo = Xo(:, [2:12]);
        Yo = Yo(:, [2:12]);
        Zo = Zo(:, [2:12]);
        X = [X Xo];
        Y = [Y Yo];
        Z = [Z Zo];
        C = [C {[0 0 0]}];
        
        [X, Y, Z] = shift_instrument(X, Y, Z, 0, 0, 2*inner);
        
        [X, Y, Z, C] = add_cyl_to_design_3(X, Y, Z, C, [-outer*0.58 0 0], [-outer*0.58 0 3.5*inner], inner , 0, [1 1], [0 0 0]);
        [X, Y, Z, C] = add_cyl_to_design_3(X, Y, Z, C, [outer*0.58 0 0], [outer*0.58 0 3.5*inner], inner , 0, [1 1], [0 0 0]);
        [X, Y, Z, C] = add_cyl_to_design_3(X, Y, Z, C, [-outer 0 inner], [outer 0 inner], inner*0.8 , 0, [1 1], [01 0.8 0]);
        
        [X, Y, Z] = shift_instrument(X, Y, Z, 0, 0, -inner);
        
        
    case {'150mm load ring' 'herc-alloy oblong links 6-10mm'}
        
        overlap_factor = 1.3;
        switch lower(instrument)
            case '150mm load ring'
                L = 0.15;
                outer = 0.05;
                inner = 0.01;
                n = 13;
            case 'herc-alloy oblong links 6-10mm'
                L = 0.1;
                outer = 0.05;
                inner = 0.01;
                n = 13;
        end
        
        L = L*overlap_factor;
        outer = outer*overlap_factor;
        inner = inner*overlap_factor;
        
        [Xo, Yo, Zo] = get_toroid_xz(outer, inner, n);
        Xo = Xo(:, [4:10]);
        Yo = Yo(:, [4:10]);
        Zo = Zo(:, [4:10]);
        X = [X Xo];
        Y = [Y Yo];
        Z = [Z Zo];
        C = [C {[01 0.8 0]}];
        
        [X, Y, Z] = shift_instrument(X, Y, Z, 0, 0, L-2*outer);
        
        [Xo, Yo, Zo] = get_toroid_xz(outer, inner, n);
        Xo = Xo(:, [10:n 1:4]);
        Yo = Yo(:, [10:n 1:4]);
        Zo = Zo(:, [10:n 1:4]);
        X = [X Xo];
        Y = [Y Yo];
        Z = [Z Zo];
        C = [C {[01 0.8 0]}];
        
        [X, Y, Z, C] = add_cyl_to_design_3(X, Y, Z, C, [-outer 0 outer], [-outer 0 L-outer], inner , 0, [1 1], [01 0.8 0]);
        [X, Y, Z, C] = add_cyl_to_design_3(X, Y, Z, C, [outer 0 outer], [outer 0 L-outer], inner , 0, [1 1], [01 0.8 0]);
        [X, Y, Z] = rotate_instrument(X, Y, Z, 90, 0, 0);
        [X, Y, Z] = shift_instrument(X, Y, Z, 0, 0, L*(1-overlap_factor)/2);
        
    case 'dbld 3/8" chain'
        
        [X, Y, Z, C] = add_instrument_to_design(X, Y, Z, C, '16mmbs', -0.09, 0, 0.1, 0, 0, 0);
        [X, Y, Z, C] = add_instrument_to_design(X, Y, Z, C, '16mmbs', 0.09, 0, 0.1, 0, 0, 0);
        
        [X, Y, Z, C] = add_instrument_to_design(X, Y, Z, C, '20cm_chain', 0.015, 0, 0, 0, 0, 37);
        [X, Y, Z, C] = add_instrument_to_design(X, Y, Z, C, '20cm_chain', -0.015, 0, 0, 0, 0, -37);
        [X, Y, Z, C] = add_instrument_to_design(X, Y, Z, C, 'linkxz', 0, 0.015, -0.015, 0, 45, 0);
    
    case 'sercel beacon'
        
        col = [1 0.8 0.2];
        [X, Y, Z, C] = add_shs_to_design_3(X, Y, Z, C, [0 0 0], [0 0 0.1], [0.1 0.05] , 0, [0.95 0.95], [1 0.8 0.2]);
        [X, Y, Z, C] = add_shs_to_design_3(X, Y, Z, C, [0 0 0.1], [0 0 0.15], [0.1 0.05] , 0, [1 0.95], [1 0.8 0.2]);
        [X, Y, Z, C] = add_shs_to_design_3(X, Y, Z, C, [0 0 0.15], [0 0 0.2], [0.1 0.05] , 0, [1 1], [0.7 0.7 0.7]);
        
        [X, Y, Z, C] = add_shs_to_design_3(X, Y, Z, C, [0 0 0.05], [0 0 0.08], [0.1 0.05] , 0, [0.95000001 0.9500001], [1 0.5 0.0]);
        [X, Y, Z, C] = add_shs_to_design_3(X, Y, Z, C, [0 0 0.16], [0 0 0.19], [0.1 0.05] , 0, [1.000001 1.000001], [1 0.5 0.0]);
        
        [X, Y, Z, C] = add_shs_to_design_3(X, Y, Z, C, [0 0 0.2], [0 0 0.23], [0.015 0.015] , 0, [1 1], [0.9 0.9 0.9]);
        [X, Y, Z, C] = add_shs_to_design_3(X, Y, Z, C, [0 0 0.23], [0 0 0.29], [0.01 0.01] , 0, [1 1], [0.5 0.5 0.5]);
        
    case {'sbe37 ctd' 'sbe37 ct'}
        
        [X, Y, Z, C] = add_cyl_to_design_3(X, Y, Z, C, [0 0 0], [0 0 0.56], 0.05 , 0, [1 1], [0.99 0.99 0.99]);
        [X, Y, Z, C] = add_shs_to_design_3(X, Y, Z, C, [-0.05 0 0.28], [-0.05 0 0.56-0.0001], [0.09 0.08] , 0, [1 1], [0.8 0.8 0.8]);
        
        [X, Y, Z, C] = add_instrument_to_design(X, Y, Z, C, 'sbe_plug', 0, 0, 0.56, 0, 0, 0);
        
        [X, Y, Z] = shift_instrument(X, Y, Z, 0, 0, -0.28);
        
    case {'sbe39 tp' 'sbe39 t'}
        
        [X, Y, Z, C] = add_cyl_to_design_3(X, Y, Z, C, [0 0 0], [0 0 0.3], 0.03 , 0, [1 1], [0.99 0.99 0.99]);
        [X, Y, Z, C] = add_cyl_to_design_3(X, Y, Z, C, [0 0 0.26], [0 0 0.31], 0.031 , 0, [1 1], [0.2 0.2 0.2]);
        
        [X, Y, Z] = shift_instrument(X, Y, Z, 0, 0, -0.15);
        
    case {'sbe39-ext'}
        
        [X, Y, Z, C] = add_instrument_to_design(X, Y, Z, C, 'sbe39 t', 0, 0, 0, 0, 0, 0);
        [X, Y, Z, C] = add_instrument_to_design(X, Y, Z, C, 'sbe_plug', 0, 0, 0.15, 0, 0, 0);
        
    case {'sbe_plug'}
        
        [X, Y, Z, C] = add_cyl_to_design_3(X, Y, Z, C, [0 0 0], [0 0 0.02], 0.011 , 0, [1 1], [0 0 0]);
        [X, Y, Z, C] = add_cyl_to_design_3(X, Y, Z, C, [0 0 0], [0 0 0.05], 0.008 , 0, [1 1], [0 0 0]);
        [X, Y, Z, C] = add_cyl_to_design_3(X, Y, Z, C, [0 0 0.05], [0 0 0.07], 0.011 , 0, [1 1], [0 0 0]);
        
    case {'vemco t'}
        
        [X, Y, Z, C] = add_cyl_to_design_3(X, Y, Z, C, [0 0 -0.2], [0 0 0.4-0.2], 0.03 , 0, [0.85 0.95], [0.2 0.2 0.2]);
    
    case {'fat arrow'}
        
        asize = 1;
        
        colour = [0.8 0.2 0.2];
        for ii = 1:2:length(varargin)
            switch lower(varargin{ii})
                case 'color'
                    switch lower(varargin{ii+1})
                        case 'blue'
                            colour = [0.2 0.2 0.8];
                        case 'green'
                            colour = [0.2 0.8 0.2];

                    end
                case 'size'
                    asize  = varargin{ii+1};
            end
        end
        [X, Y, Z, C] = add_cyl_to_design_3(X, Y, Z, C, [0 0 0]*asize , [0.3 0 0]*asize , 0.03*asize  , 1*asize , [1 1], colour);
        [X, Y, Z, C] = add_cyl_to_design_4(X, Y, Z, C, [0.3 0 0]*asize , [0.5 0 0]*asize , 0.1*asize  , 0*asize , [1 1], [1 0], colour);
            
    otherwise
        
         fprintf(1, 'No drawing data for %s', instrument)
         error
        
end

lpp=inf;
lpn=inf;

if inv
    rotmat = rotation_matrix(0, 0, 180);
    for ii = 1:numel(Z)
        lpp = min([lpp min(Z{ii}(:))]);
        n = size(X{ii}, 1);
        vector = surf2vector(X{ii}, Y{ii}, Z{ii});
        vector = rotmat*vector;
        [X{ii}, Y{ii}, Z{ii}] = vector2surf(vector, n);
        lpn = min([lpn min(Z{ii}(:))]);
    end
    for ii = 1:numel(Z)
        Z{ii} = Z{ii} + lpp-lpn;
    end
end


end

function [x, y, z] = cyl_a2b_reprise(a, b, d)

[cylX, cylY, cylZ] = cylinder;
[sphX,sphY,sphZ] = sphere;

v = b - a;
l = norm(v);

heading = -180*atan2(v(2), v(1))/pi; % - dx/dy
pitch = 0; % about x
roll = 180*atan2(sqrt(v(1)^2+v(2)^2), v(3))/pi; % about y - dx/dz
% heading = 0;
rotmat = rotation_matrix(heading, pitch, roll, [1 2 3]);
vector = surf2vector(cylX*d, cylY*d, cylZ*l);
vector = rotmat*vector;
[x, y, z] = vector2surf(vector, 2);
x = {x + a(1)};
y = {y + a(2)};
z = {z + a(3)};

x{end+1} =  sphX*d+a(1);
z{end+1} =  sphZ*d+a(3);
y{end+1} =  sphY*d+a(2);

x{end+1} =  sphX*d+b(1);
z{end+1} =  sphZ*d+b(3);
y{end+1} =  sphY*d+b(2);

end

function [Xn, Yn, Zn] = cyl_a2b_reprise_2(a, b, d, zf)

[cylX, cylY, cylZ] = cylinder;
[sphX,sphY,sphZ] = sphere;

v = b - a;
l = norm(v);

heading = -180*atan2(v(2), v(1))/pi; % - dx/dy
pitch = 0; % about x
roll = 180*atan2(sqrt(v(1)^2+v(2)^2), v(3))/pi; % about y - dx/dz
% heading = 0;
rotmat = rotation_matrix(heading, pitch, roll, [1 2 3]);

vector = surf2vector(cylX*d, cylY*d, cylZ*l);
vector = rotmat*vector;
[x, y, z] = vector2surf(vector, 2);
Xn = {x + a(1)};
Yn = {y + a(2)};
Zn = {z + a(3)};

n = size(sphX, 1);
vector = surf2vector(sphX*d, sphY*d, sphZ*d*zf);
vector = rotmat*vector;
[x, y, z] = vector2surf(vector, n);
Xn = cat(2, Xn, {x + a(1)});
Yn = cat(2, Yn, {y + a(2)});
Zn = cat(2, Zn, {z + a(3)});

vector = surf2vector(sphX*d, sphY*d, sphZ*d*zf+l);
vector = rotmat*vector;
[x, y, z] = vector2surf(vector, n);
Xn = cat(2, Xn, {x + a(1)});
Yn = cat(2, Yn, {y + a(2)});
Zn = cat(2, Zn, {z + a(3)});

end

function [Xn, Yn, Zn] = shs_a2b_reprise_3(a, b, d, zf, tf)

%[cylX, cylY, cylZ] = cylinder;
%[sphX,sphY,sphZ] = sphere;
shsX = [-1 -1 1 1 -1 ;
        -1 -1 1 1 -1];
shsY = [-1 1 1 -1 -1 ;
        -1 1 1 -1 -1];
shsZ = [1 1 1 1 1 ;
        0 0 0 0 0];

v = b - a;
l = norm(v);

heading = -180*atan2(v(2), v(1))/pi; % - dx/dy
pitch = 0; % about x
roll = 180*atan2(sqrt(v(1)^2+v(2)^2), v(3))/pi; % about y - dx/dz
% heading = 0;
rotmat = rotation_matrix(heading, pitch, roll, [1 2 3]);

% vector = surf2vector(shsX*d(1), shsY*d(2), shsZ*l);

shsX(1, :) = shsX(1, :)*(d(1)/2)*tf(1);
shsX(2, :) = shsX(2, :)*(d(1)/2)*tf(2);
shsY(1, :) = shsY(1, :)*(d(2)/2)*tf(1);
shsY(2, :) = shsY(2, :)*(d(2)/2)*tf(2);
shsZ = shsZ*l;

vector = surf2vector(shsX, shsY, shsZ);
vector = rotmat*vector;
[x, y, z] = vector2surf(vector, 2);
if exist('tf', 'var')
    x(1, :) = x(1, :) * tf(1);
    x(2, :) = x(2, :) * tf(2);
    y(1, :) = y(1, :) * tf(1);
    y(2, :) = y(2, :) * tf(2);
end
x = x + a(1);
y = y + a(2);
z = z + a(3);

Xn = {x};
Yn = {y};
Zn = {z};

Xn{2} = x(1, :);
Yn{2} = y(1, :);
Zn{2} = z(1, :);

Xn{3} = x(2, :);
Yn{3} = y(2, :);
Zn{3} = z(2, :);

end

function [Xn, Yn, Zn] = cyl_a2b_reprise_3(a, b, d, zf, tf, tfx)

[cylX, cylY, cylZ] = cylinder;
[sphX,sphY,sphZ] = sphere;

if numel(zf) == 1
    zf = zf*[1 1];
end

v = b - a;
l = norm(v);

heading = -180*atan2(v(2), v(1))/pi; % - dx/dy
pitch = 0; % about x
roll = 180*atan2(sqrt(v(1)^2+v(2)^2), v(3))/pi; % about y - dx/dz
% heading = 0;
rotmat = rotation_matrix(heading, pitch, roll, [1 2 3]);

vector = surf2vector(cylX*d, cylY*d, cylZ*l);
vector = rotmat*vector;
[x, y, z] = vector2surf(vector, 2);
if exist('tf', 'var')
    x(1, :) = x(1, :) * tf(1);
    x(2, :) = x(2, :) * tf(2);
    y(1, :) = y(1, :) * tf(1);
    y(2, :) = y(2, :) * tf(2);
end
if exist('tfx', 'var')
    z(1, :) = z(1, :) * tfx(1);
    z(2, :) = z(2, :) * tfx(2);
    y(1, :) = y(1, :) * tfx(1);
    y(2, :) = y(2, :) * tfx(2);
end
Xn = {x + a(1)};
Yn = {y + a(2)};
Zn = {z + a(3)};

n = size(sphX, 1);
vector = surf2vector(sphX*d, sphY*d, sphZ*d*zf(1));
vector = rotmat*vector;
[x, y, z] = vector2surf(vector, n);
if exist('tf', 'var')
    x(:, :) = x(:, :) * tf(1);
    y(:, :) = y(:, :) * tf(1);
end
if exist('tfx', 'var')
    z(:, :) = z(:, :) * tfx(1);
    y(:, :) = y(:, :) * tfx(1);
end
Xn = cat(2, Xn, {x + a(1)});
Yn = cat(2, Yn, {y + a(2)});
Zn = cat(2, Zn, {z + a(3)});

vector = surf2vector(sphX*d, sphY*d, sphZ*d*zf(2)+l);
vector = rotmat*vector;
[x, y, z] = vector2surf(vector, n);
if exist('tf', 'var')
    x(:, :) = x(:, :) * tf(2);
    y(:, :) = y(:, :) * tf(2);
end
if exist('tfx', 'var')
    z(:, :) = z(:, :) * tfx(2);
    y(:, :) = y(:, :) * tfx(2);
end
Xn = cat(2, Xn, {x + a(1)});
Yn = cat(2, Yn, {y + a(2)});
Zn = cat(2, Zn, {z + a(3)});

end

function [X, Y, Z, C] = add_inverted_instrument_to_design(X, Y, Z, C, name, xo, yo, zo, r1, r2, r3)

    [x, y, z, c] = draw_instrument(name, true);
    [x, y, z] = rotate_instrument(x, y, z, r1, r2, r3);
    [x, y, z] = shift_instrument(x, y, z, xo, yo, zo);
    [X, Y, Z, C] = add_XYZC_to_design(X, Y, Z, C, x, y, z, c);

end

function [X, Y, Z, C] = add_instrument_to_design(X, Y, Z, C, name, xo, yo, zo, r1, r2, r3)

    [x, y, z, c] = draw_instrument(name, false);
    [x, y, z] = rotate_instrument(x, y, z, r1, r2, r3);
    [x, y, z] = shift_instrument(x, y, z, xo, yo, zo);
    [X, Y, Z, C] = add_XYZC_to_design(X, Y, Z, C, x, y, z, c);

end

function [X, Y, Z, C] = add_XYZC_to_design(X, Y, Z, C, x, y, z, c)

X = cat(2, X, x);
Y = cat(2, Y, y);
Z = cat(2, Z, z);
if iscell(c)
    if ~isequal(size(c), size(z))
        c = repmat({c}, 1, numel(x));
    end
else
    c = repmat({c}, 1, numel(x));
end
C = cat(2, C, c);  

end

function [X, Y, Z, C] = add_cyl_to_design(X, Y, Z, C, a, b, d ,c)

[x, y, z] = cyl_a2b_reprise(a, b, d);
[X, Y, Z, C] = add_XYZC_to_design(X, Y, Z, C, x, y, z, c);

end

function [X, Y, Z, C] = add_cyl_to_design_2(X, Y, Z, C, a, b, d, zf ,c)

[x, y, z] = cyl_a2b_reprise_2(a, b, d, zf);
[X, Y, Z, C] = add_XYZC_to_design(X, Y, Z, C, x, y, z, c);

end

function [X, Y, Z, C] = add_cyl_to_design_3(X, Y, Z, C, a, b, d, zf, tf, c)

if any(c>1)
    error('Bad colour')
end
[x, y, z] = cyl_a2b_reprise_3(a, b, d, zf, tf);
[X, Y, Z, C] = add_XYZC_to_design(X, Y, Z, C, x, y, z, c);

end

function [X, Y, Z, C] = add_cyl_to_design_4(X, Y, Z, C, a, b, d, zf, tf, ftx,c)

[x, y, z] = cyl_a2b_reprise_3(a, b, d, zf, tf, ftx);
[X, Y, Z, C] = add_XYZC_to_design(X, Y, Z, C, x, y, z, c);

end

function [X, Y, Z, C] = add_shs_to_design_3(X, Y, Z, C, a, b, d, zf, tf ,c)

[x, y, z] = shs_a2b_reprise_3(a, b, d, zf, tf);
[X, Y, Z, C] = add_XYZC_to_design(X, Y, Z, C, x, y, z, c);

end

function [x, y, z] = get_toroid_xz(r1, r2, n)

theta  = linspace(-2*pi, 0, n)   ; % Poloidal angle
phi    = linspace(-pi/2, 3.*pi/2, n) ; % Toroidal angle
[t, p] = meshgrid(phi, theta);
x = (r1 + r2.*cos(p)) .* cos(t);
z = (r1 + r2.*cos(p)) .* sin(t);
y = r2.*sin(p);

z = z + r1;
end

function [x, y, z] = get_toroid_yz(r1, r2, n)

theta  = linspace(-pi, pi, n)   ; % Poloidal angle
phi    = linspace(0., 2.*pi, n) ; % Toroidal angle
[t, p] = meshgrid(phi, theta);
y = (r1 + r2.*cos(p)) .* cos(t);
z = (r1 + r2.*cos(p)) .* sin(t);
x = r2.*sin(p);

z = z + r1;
end
