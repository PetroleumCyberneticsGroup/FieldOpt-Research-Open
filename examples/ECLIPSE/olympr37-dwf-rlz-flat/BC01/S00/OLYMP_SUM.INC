-- ---------------------------------------------------------
RUNSUM     -- additional table in *.PRT file
--SEPARATE   -- write a seperate *.RSM file
NARROW     -- increases the number of columns per page
RPTONLY    -- summary data only produced at report times

--RPTSMRY
--1 /

SUMTHIN
50 /

-- PAGE 1 --------------------------------------------------
TIMESTEP--timestep length                                  1
FPR     -- field pressure                                  2
FPRH    -- field pressure by hydrocarbon volume            3
FRPV    -- pore volume at reservoir conditions              4
FOE     -- oil recovery                                    5
FOIP    -- field oil in place                              6
FWIP    -- field water in place                            7
FGIP    -- field gas in place (SPACING) .................. 8
ELAPSED -- elapsed time in seconds                         9
MEMORYTS-- max current mem usage across procs os-reported 10
TCPU    -- current cpu usage in secs                      11

-- PAGE 2 --------------------------------------------------
-- FMWPR   -- # of prods currently flowing
-- FMWPG   -- # of prods controlled on group control
-- FMWPO   -- # of prods controlled by own oil rate limit
-- FMWPP   -- # of prods controlled on pressure control
-- FMWIN   -- # of injs currently flowing

-- PAGE X --------------------------------------------------
SOFR    -- segment oil flow rate
PRODX2 /
/
FGIP    -- SPACING ......................................
FGIP    -- SPACING ......................................
FGIP    -- SPACING ......................................
FGIP    -- SPACING ......................................
FGIP    -- SPACING ......................................
FGIP    -- SPACING ......................................
FGIP    -- SPACING ......................................
FGIP    -- SPACING ......................................

-- PAGE X --------------------------------------------------
SWFR    -- segment water flow rate
PRODX2 /
/
FGIP    -- SPACING ......................................
FGIP    -- SPACING ......................................
FGIP    -- SPACING ......................................
FGIP    -- SPACING ......................................
FGIP    -- SPACING ......................................
FGIP    -- SPACING ......................................
FGIP    -- SPACING ......................................
FGIP    -- SPACING ......................................

-- PAGE X --------------------------------------------------
SGFR    -- segment water flow rate
PRODX2 /
/
FGIP    -- SPACING ......................................
FGIP    -- SPACING ......................................
FGIP    -- SPACING ......................................
FGIP    -- SPACING ......................................
FGIP    -- SPACING ......................................
FGIP    -- SPACING ......................................
FGIP    -- SPACING ......................................
FGIP    -- SPACING ......................................

-- PAGE X --------------------------------------------------
SPR     -- segment pressure
PRODX2 /
/
FGIP    -- SPACING ......................................
FGIP    -- SPACING ......................................
FGIP    -- SPACING ......................................
FGIP    -- SPACING ......................................
FGIP    -- SPACING ......................................
FGIP    -- SPACING ......................................
FGIP    -- SPACING ......................................
FGIP    -- SPACING ......................................

-- PAGE X --------------------------------------------------
SPRD    -- segment pressure drop
PRODX2 /
/
FGIP    -- SPACING ......................................
FGIP    -- SPACING ......................................
FGIP    -- SPACING ......................................
FGIP    -- SPACING ......................................
FGIP    -- SPACING ......................................
FGIP    -- SPACING ......................................
FGIP    -- SPACING ......................................
FGIP    -- SPACING ......................................

-- PAGE X --------------------------------------------------
SWCT    -- segment water cut
PRODX2 /
/
FGIP    -- SPACING ......................................
FGIP    -- SPACING ......................................
FGIP    -- SPACING ......................................
FGIP    -- SPACING ......................................
FGIP    -- SPACING ......................................
FGIP    -- SPACING ......................................
FGIP    -- SPACING ......................................
FGIP    -- SPACING ......................................

-- PAGE X --------------------------------------------------
-- segment cross sectional area (e300)
SCSA
PRODX2 /
/
FGIP    -- SPACING ......................................
FGIP    -- SPACING ......................................
FGIP    -- SPACING ......................................
FGIP    -- SPACING ......................................
FGIP    -- SPACING ......................................
FGIP    -- SPACING ......................................
FGIP    -- SPACING ......................................
FGIP    -- SPACING ......................................

-- PAGE X --------------------------------------------------
-- COPR    -- oil flow rate in lumped completion
-- PRODX2 /
-- /
-- FGIP    -- SPACING ......................................
-- FGIP    -- SPACING ......................................

--CWPR    -- water flow rate in lumped completion
--PRODX2 /
--/
--FGIP    -- SPACING ......................................
--FGIP    -- SPACING ......................................

--CPR     -- connection pressure
--PRODX2 /
--/
--FGIP    -- SPACING ......................................
--FGIP    -- SPACING ......................................

--CPI     -- productivity index of well's preferred phase
--PRODX2 /
--/
--FGIP    -- SPACING ......................................
--FGIP    -- SPACING ......................................

--CTFAC   -- connection transmissibility factor
--PRODX2 /
--/
--FGIP    -- SPACING ......................................
--FGIP    -- SPACING ......................................

-- PAGE X --------------------------------------------------
-- BHP | STATUS | CONTROL

WBHP
'PROD*' /

WBHP
'INJ*' /

-- -------
WSTAT
'INJ*' /

WMCTL
'INJ*' /

WSTAT
'PROD*' /

WMCTL
'PROD*' /

FGIP    -- SPACING ......................................
FGIP    -- SPACING ......................................
FGIP    -- SPACING ......................................

-- PAGE X --------------------------------------------------
-- FIELD TOTALS

FOPT

FWPT

FGPT

FLPT

FWIT

FGIT

FWCT

FGOR

FGIP    -- SPACING ......................................
FGIP    -- SPACING ......................................
FGIP    -- SPACING ......................................
FGIP    -- SPACING ......................................

-- PAGE X --------------------------------------------------
-- FIELD RATES

FOPR

FWPR

FGPR

FLPR

FWIR

FGIR

FGIP    -- SPACING ......................................
FGIP    -- SPACING ......................................
FGIP    -- SPACING ......................................
FGIP    -- SPACING ......................................
FGIP    -- SPACING ......................................
FGIP    -- SPACING ......................................

-- PAGE X --------------------------------------------------
-- WELL TOTALS

WOPT
'PROD*' /

WWPT
'PROD*' /

WGPT
'PROD*' /

WLPT
'PROD*' /

WWIT
'PROD*'
'INJ*' /

WGIT
'PROD*'
'INJ*' /

-- RATIOS
WWCT
'PROD*' /

WGOR
'PROD*' /

-- PAGE X --------------------------------------------------
-- WELL RATES: PRODX

WOPR
'PROD*' /

WWPR
'PROD*' /

WGPR
'PROD*' /

WLPR
'PROD*' /

WWIR
'PROD*'
'INJ*' /

WGIR
'PROD*'
'INJ*' /

FMCTP

FGIP    -- SPACING ......................................
