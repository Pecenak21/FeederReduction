function pen = getPenLevel(c)
pen = sum([c.pvsystem.Pmpp])/ sum([c.load.kw])*100;
end