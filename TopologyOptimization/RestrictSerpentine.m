function flag=RestrictSerpentine(NODE,BARS)
tol = 1e-3;

flag = rCircle([2 10.9],10.1-tol,NODE,BARS) | ...
       rCircle([6 -8.9],10.1-tol,NODE,BARS);