function c=muted_colors(ind)
Palette= [    [51 34 136];... % indigo
                    [136 204 238];... % cyan
                    [68 170 153];... % teal
                    [17 119 51];... % green
                    [153 153 51];... % olive
                    [221 204 119];... % sand
                    [204 102 119];... % rose
                    [136 34 85];... % wine
                    [170 68 153];... % purple
    ]./256;
c=Palette(ind,:);
end