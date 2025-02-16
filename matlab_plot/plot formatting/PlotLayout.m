function P=PlotLayout(fhs,x,y)

width=sum(x);
height=sum(y);

rows=2:2:length(y);
cols=2:2:length(x);

nRows=length(rows);
nCols=length(cols);

P=cell(nRows,nCols);

for rr=1:nRows
    for cc=1:nCols
        fromLeft=sum(x(1:cols(cc)-1))/width;
        fromBott=sum(y(1:rows(rr)-1))/height;
        w=x(cols(cc))/width;
        h=y(rows(rr))/height;
        P{rr,cc}=[fromLeft fromBott w h ];
    end
end

P=P(end:-1:1,:);

Pix_SS = get(0,'screensize');

for ff=1:length(fhs)
    fh=fhs(ff);
    
    set(fh,'PaperUnits','centimeters')
    set(fh,'PaperSize',[width height])
    set(fh,'PaperPosition',[0 0 width height])

    units=get(fh,'units');
    
    set(fh,'units','centimeters')
    pos=get(fh,'position');pos1 = pos;
    pos(3:4)=[width height];pos(2) = pos1(2) + pos1(4) - pos(4);
    set(fh,'position',pos);
    %% Screensize
    set(fh,'units','pixels');
    pos=get(fh,'position');poso = pos;
    if Pix_SS(3)<(sum(pos([1 3])))+70
        pos(1)=pos(1)-(sum(pos([1 3])))-70+Pix_SS(3);
    end
    if Pix_SS(4)<(sum(pos([2 4])))+50
        pos(2)=pos(2)-(sum(pos([2 4])))-70+Pix_SS(4);
    end
    if sum(abs(poso(3:4) - pos(3:4)))
        set(fh,'position',pos);
    end
    %%
    set(fh,'units',units);
end


end