ACTIONX
    ICVS_T0  1  0  /
    DAY  >=  0  /
/

WSEGVALV
--WName  SegNo  Cv   Ac         L   D       R   Pac         Stat  Amax     -
PRODX2   2      0.7  3.9265E-3  1*  0.1016  1*  8.10732E-3  OPEN  7.853E-3 /
PRODX2   3      0.7  3.9265E-3  1*  0.1016  1*  8.10732E-3  OPEN  7.853E-3 /
PRODX2   4      0.7  3.9265E-3  1*  0.1016  1*  8.10732E-3  OPEN  7.853E-3 /
PRODX2   5      0.7  3.9265E-3  1*  0.1016  1*  8.10732E-3  OPEN  7.853E-3 /
PRODX2   6      0.7  3.9265E-3  1*  0.1016  1*  8.10732E-3  OPEN  7.853E-3 /
PRODX2   7      0.7  3.9265E-3  1*  0.1016  1*  8.10732E-3  OPEN  7.853E-3 /
PRODX2   8      0.7  3.9265E-3  1*  0.1016  1*  8.10732E-3  OPEN  7.853E-3 /
PRODX2   9      0.7  3.9265E-3  1*  0.1016  1*  8.10732E-3  OPEN  7.853E-3 /
/

ENDACTIO

--item 1:  Well name,
--item 2:  Segment number to contain the valve
--item 3:  Dimensionless flow coefficient for the valve, Cv
--item 4:  Cross-section area for flow in the constriction, Ac

--item 5:  Additional length of pipe for the friction pressure drop, L
--item 6:  Pipe diameter, D (for frictional pressure drop calculation)
--item 7:  Absolute roughness of pipe wall.

--itme 8:  Pipe cross-section area, Ap (for frictional pressure drop calculation++)
--itme 9:  Flag to indicate whether device is open or shut
--itme 10: Maximum cross-sectional area for flow in the constriction, Amax