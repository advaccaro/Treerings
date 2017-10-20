for i =1:512
    plot(tcurve(i).trdat(10:end),tcurve(i).stdat(1:end-9),'*');
end