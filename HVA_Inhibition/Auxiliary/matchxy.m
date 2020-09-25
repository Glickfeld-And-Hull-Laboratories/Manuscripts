function matchxy(minType,dotLine)
tempMax = max([xlim,ylim]);

if strcmp(minType,'min')
    tempMin = min([xlim,ylim]);
    xlim([tempMin tempMax]); ylim(xlim);
else
    tempMin = 0;
    xlim([tempMin tempMax]); ylim(xlim);
end

if dotLine
    hold on
    plot([tempMin tempMax],[tempMin tempMax],'k--');
end

