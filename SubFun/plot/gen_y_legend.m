function Y_legend = gen_y_legend(Y_legend, DATA, nb_subjects)
%  generate a legend for the Y axis when displaying examples of models for
%  Linear mixed models and F test

Y_legend{end+1} = [DATA.Name ' ; sub-01 ; CST' ]; %#ok<*AGROW>
for i=1:nb_subjects-2
    Y_legend{end+1} = '...';
end
Y_legend{end+1} = [DATA.Name ' ; sub-11 ; CST' ];
Y_legend{end+1} = [DATA.Name ' ; sub-01 ; LIN' ];
for i=1:nb_subjects-2
    Y_legend{end+1} = '...';
end
Y_legend{end+1} = [DATA.Name ' ; sub-11 ; LIN' ];
end