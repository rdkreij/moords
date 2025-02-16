function [ output_args ] = actually_linkaxes( axes, aa )

    error('this doesn''t work - endless loops')
    
    addListeners

    function addListeners
        % Create listeners. We store the array of listeners in the axes to make
        % sure that they have the same life-span as the axes they are listening to.
        
        ax = axes(aa);
        
        axh = handle( ax );
        listeners = [
            addlistener( axh, findprop( axh, 'XLim' ), 'PostSet', @execute_this )
            addlistener( axh, findprop( axh, 'YLim' ), 'PostSet', @execute_this )
            ];
        setappdata( ax, 'actually_linkaxes_listeners', listeners );
        setappdata( ax, 'actually_linkaxes_linkedaxes', axes );
        
    end % addListeners

    function removeListeners
        % Rempove any property listeners whilst we are fiddling with the axes
        if isappdata( ax, 'actually_linkaxes_listeners' )
            delete( getappdata( ax, 'actually_linkaxes_listeners' ) );
            rmappdata( ax, 'actually_linkaxes_listeners' );
        end
    end % removeListeners

    function execute_this(varargin)
        
        removeListeners
        
        varargin{1};
        ax = axes(aa);
        
        xl = get(ax, 'xlim');
        yl = get(ax, 'ylim');
        
        set(axes, 'xlim', xl)
        set(axes, 'ylim', yl)
       
        addListeners
        
    end

end

