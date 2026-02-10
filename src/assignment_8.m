% Adapted from https://se.mathworks.com/help/pde/ug/constant-boundary-condition-specifications.html
epsilon = 0.1;

C1 = [1,0,0,epsilon]';
C2 = [1,0,0,1]';

C1 = [C1;zeros(length(C2)-length(C1),1)];
geom = [C2,C1];

ns = (char('C2','C1'))';

sf = 'C2 - C1';

g = decsg(geom,sf,ns);

model = createpde;

geometryFromEdges(model,g);
pdegplot(model,EdgeLabels="on")
xlim([-1.1 1.1])
axis equal