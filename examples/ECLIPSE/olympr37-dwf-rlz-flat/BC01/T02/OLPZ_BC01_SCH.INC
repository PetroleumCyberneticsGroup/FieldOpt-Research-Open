--OLYMPUS_37_SCH_OPT.INC

-- _________________________________________________________
INCLUDE
'OLPZ_XICD_WLS.INC' /

-- _________________________________________________________
WTEST
'PROD*' 1  P 1* /
'INJD*' 1  P 1* /

/

GRUPTREE
'P' FIELD /
'I' FIELD /
/

WCONPROD
PRODX2  1* BHP 5* 155.0000 /
/

--ORIG: 235
WCONINJE
INJD-15 WATER 1* BHP 1* 1* 205.0000 /
INJD-16 WATER 1* BHP 1* 1* 205.0000 /
/
-- _________________________________________________________
INCLUDE
'fo_edits.INC'/

-- _________________________________________________________
TSTEP

--6 years
200*10.95 /

