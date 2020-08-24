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

-- PAGE X --------------------------------------------------
COPR    -- oil flow rate in lumped completion
PRODX2 /
/
FGIP    -- SPACING ......................................
FGIP    -- SPACING ......................................
FGIP    -- SPACING ......................................
FGIP    -- SPACING ......................................
FGIP    -- SPACING ......................................

CWPR    -- water flow rate in lumped completion
PRODX2 /
/
FGIP    -- SPACING ......................................
FGIP    -- SPACING ......................................
FGIP    -- SPACING ......................................
FGIP    -- SPACING ......................................
FGIP    -- SPACING ......................................

CPR     -- connection pressure
PRODX2 /
/
FGIP    -- SPACING ......................................
FGIP    -- SPACING ......................................
FGIP    -- SPACING ......................................
FGIP    -- SPACING ......................................
FGIP    -- SPACING ......................................

CPI     -- productivity index of well's preferred phase
PRODX2 /
/

CTFAC   -- connection transmissibility factor
PRODX2 /
/

-- PAGE 2 --------------------------------------------------

WBHP
'INJ*' /

WSTAT
'INJ*' /

WMCTL
'INJ*' /

-- PAGE 3+4 ------------------------------------------------

WBHP
'PROD*' /

WSTAT
'PROD*' /

WMCTL
'PROD*' /

-- PAGE 3 --------------------------------------------------



-- STATES (WELLS)

WSTAT
'PROD*' /




FMCTP

-- RATES (WELLS)

WOPR
'PROD*' /

WWPR
'PROD*' /

WLPR
'PROD*' /

WGPR
'PROD*' /

WWIR
'INJ*' /

-- RATES (FIELD)

FOPR

FWPR

FLPR

FGPR

FWIR

-- TOTALS (WELLS)

WOPT
'PROD*' /

WWPT
'PROD*' /

WLPT
'PROD*' /

WGPT
'PROD*' /

WWIT
'INJ*' /

-- TOTALS (FIELD)

FOPT

FWPT

FLPT

FGPT

FWIT


-- RATIOS
WWCT
'PROD*' /

FWCT

WGOR
'PROD*' /

FGOR

