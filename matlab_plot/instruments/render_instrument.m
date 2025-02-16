function render_instrument( X, Y, Z, C, ah )

if nargin == 4
    ah = gca;
end

ec = 'none';
hold on;
for cc = 1:numel(X)
    if isvector(X{cc})
        ph = fill3(X{cc}, Y{cc}, Z{cc}, 'k' ,'facecolor', C{cc}, 'edgecolor', ec, 'parent', ah);
    else
        try
            ph = surf(ah, X{cc}, Y{cc}, Z{cc} ,'facecolor', C{cc}, 'edgecolor', ec);
        catch
            hi=1
        end
    end
end
   
axis equal
light               % create a light
lighting gouraud    % prefered method for lighting curved surfaces

end

