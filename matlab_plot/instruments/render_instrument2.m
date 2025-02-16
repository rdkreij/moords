function render_instrument2( X, Y, Z, C, scale, top, ec, fc, axish, name, list)

one = false;
if iscell(X)
    for cc = 1:numel(X)
        if isvector(X{cc})
            ph = fill3(X{cc}*scale, Y{cc}*scale, Z{cc}*scale+top, 'k', 'facecolor', C{cc}, 'edgecolor', ec, 'parent', axish);
            done = true;
        else
            ph = surf(axish, X{cc}*scale, Y{cc}*scale, Z{cc}*scale+top, 'facecolor', C{cc}, 'edgecolor', ec);
            done = true;
        end
    end
else
    ph = surf(axish, X ,Y ,Z*scale+top , 'facecolor', fc, 'edgecolor', ec);
    done = true;
%     fill3(X(2,:),Y(2,:),Z(2,:), 'k', 'facecolor', fc)
end

axis equal
% light               % create a light
% lighting gouraud    % prefered method for lighting curved surfaces

end

