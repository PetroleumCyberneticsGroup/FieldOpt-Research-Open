
OPTDIMS
100 100 /

OPTFUNC
'FOPT'   FIELD   377.3887e0  0.0  0 200 /
'FWPT'   FIELD   -37.7388e0  0.0  0 200 /
'FWIT'   FIELD   -12.5796e0  0.0  0 200 /
/

--OPTCONS
--WWCT PROD1 < 1.0 /
--/

OPTPARS
--1   2             3    4          5  6   7
--WBHP  PROD1       70   150        1  10  /
'SCSA'  'PROD1:2'   0.0  7.853E-3   0  1   BHP  /
'SCSA'  'PROD1:3'   0.0  7.853E-3   0  1   BHP  /
'SCSA'  'PROD1:4'   0.0  7.853E-3   0  1   BHP  /
'SCSA'  'PROD1:5'   0.0  7.853E-3   0  1   BHP  /
'SCSA'  'PROD1:6'   0.0  7.853E-3   0  1   BHP  /
'SCSA'  'PROD1:7'   0.0  7.853E-3   0  1   BHP  /
'SCSA'  'PROD1:8'   0.0  7.853E-3   0  1   BHP  /
'SCSA'  'PROD1:9'   0.0  7.853E-3   0  1   BHP  /
/

OPTOPTS
2 2 /
-- item 1: Steepest descent=1, Conjugate gradient=2
-- item 2: parameter dependency to OPTPARS; 1=yes, 2 (def):no

OPTTUNE
--1     2     3   4        5        6        7        8        9       10
0.05  0.05  10  1.0e-12  1.0e-12  1.0e-12  1.0e-6   1.0e-9   1.0e-6  1.0e+0 / A02: (EXP0)
--0.05  0.05  10  1.0e-6   1.0e-6   1.0e-6   1.0e-3   1.0e-6   1.0e-3  1.0e+0 / A00: (DEFAULT)
--0.05  0.05  10  1.0e-12  1.0e-12  1.0e-12  1.0e-6   1.0e-6   1.0e-3  1.0e+0 / A01: (WORK)
--0.05  0.05  10  1.0e-12  1.0e-12  1.0e-12  1.0e-6   1.0e-9   1.0e-6  1.0e+0 / A02: (EXP0)
--0.05  0.05  10  1.0e-12  1.0e-12  1.0e-12  1.0e-6   1.0e-12  1.0e-9  1.0e+0 / A03: (EXP1)
--0.05  0.05  10  1.0e-12  1.0e-12  1.0e-12  1.0e-6   1.0e-12  1.0e-9  1.0e+1 / A04: (EXP2)

-- item 1:  initial constraint tolerance (def.: 0.05)
-- item 2:  minimun constraint tolerance (def.: 0.05)
-- item 3:  # opt.iters for progression (def.: 10)

-- item 4:  obj.f. convergence tolerance (def.: 1.0e-6)
-- item 5:  minimun obj.f. RMS gradient (def.: 1.0e-6)

-- item 6:  param conv.tol b/e opt.iters (def.: 1.0e-6)
-- item 7:  param boundary rel.tol (def.: 1.0e-3)

-- item 8:  line search alpha conv.tol (def.: 1.0e-6)
-- item 9:  line search alpha minimum (def.: 1.0e-3)
-- item 10: line search alpha maximum (def.: 1.0)
