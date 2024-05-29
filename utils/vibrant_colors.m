function c=vibrant_colors(ind)
Palette=[    [0 119 187];... % blue
                    [51 187 238];... % cyan
                    [ 0 153 136];... % teal
                    [238 119 51];... % orange
                    [204 51 17];... % red
                    [238 51 119];... % magenta
                    [187 187 187];... % grey
]./256;
c=Palette(ind,:);
end