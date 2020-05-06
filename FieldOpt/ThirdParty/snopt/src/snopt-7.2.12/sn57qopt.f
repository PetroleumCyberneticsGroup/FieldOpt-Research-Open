!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!     File  sn57qopt.f
!
!     s5solv
!     s5dflt   s5Map   s5sLP    s5sQP   s5sQN   s5Stat
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s5solv
     &   ( iExit, Solver, Start,
     &     Hprod, Hprod1, QPlog, Hcalls, gotR,
     &     m, n, nb, nnH0, nnH, nName, ngQP, ngObj0, ngObj,
     &     iObj, ObjAdd, ObjQP, ObjTru, nInf, sInf,
     &     ne , nlocA, locA, indA, Acol,
     &     bl, bu, gObj, Names,
     &     nrhs0, nrhs, rhs, lenx0, nx0, x0,
     &     hEtype, hs, x, pi, rc, nS,
     &     cu, lencu, iu, leniu, ru, lenru,
     &     cw, lencw, iw, leniw, rw, lenrw )

      implicit
     &     none
      external
     &     Hprod, Hprod1, QPlog
      logical
     &     gotR
      integer
     &     Hcalls, iExit, iObj, lencu, leniu, lenru, lencw, leniw,
     &     lenrw, m, n, nb, ne, ngObj0, ngObj, ngQP, nInf, nName, nnH0,
     &     nnH, nlocA, nrhs0, nrhs, lenx0, nx0, nS, Start, locA(nlocA),
     &     indA(ne), hEtype(nb), hs(nb), iu(leniu), iw(leniw)
      double precision
     &     ObjAdd, ObjQP, ObjTru, sInf, Acol(ne), rhs(nrhs0),
     &     bl(nb), bu(nb), gObj(ngObj0), x0(lenx0), x(nb), pi(m),
     &     rc(nb), ru(lenru), rw(lenrw)
      character
     &     Solver*6, Names(nName)*8, cu(lencu)*8, cw(lencw)*8

!     ==================================================================
!     s5solv solves the current problem.
!
!     On entry,
!     the SPECS file has been read,
!     all data items have been loaded (including Acol, indA, locA, ...),
!     and workspace has been allocated within cw, iw and rw.
!     Start = lvlSrt from s3argQ.
!
!     On exit,
!     iExit  =  0 if an optimal solution was found,
!            =  1 if the problem was infeasible,
!            =  2 if the problem was unbounded,
!            =  3 if the Iteration limit was exceeded,
!           ge  4 if iterations were terminated by some other
!                 error condition (see the SQOPT user's guide).
!
!     01 Oct 1994: First version of s5solv.
!     06 Aug 1996: Min Sum option added.
!     14 Jul 1997: Thread-safe version.
!     02 Aug 2003: snEXIT and snPRNT adopted.
!     08 Mar 2004: Hot starts implemented.
!     16 May 2006: Explicit target itQP added.
!     16 Jun 2008: Hcalls added as argument.
!     18 Jun 2008: Added space for iy2, pBS and rg2.
!     ==================================================================
      character
     &     str*132, str2*132, mProb*8, istate(3)*4
      logical
     &     infsbl, needB, useQP
      integer
     &     cgItn, cgItns, eigH, gotFac, gotHes, inform, inewB,
     &     itn, itnlim, itQP, itQPmax, j, k, lAscal, lblBS,lbuBS,
     &     lblSav, lbuSav, lEmode, lenR, lgBS, lgQP, lHdx, lhEsta,
     &     lhfeas, linesL, linesS, liy, liy1, liy2, lkBS, lkx, lpBS,
     &     lR, lrg, lrg2, lsSave, lvlInf, lvlScl, lxBS,
     &     ly, ly1, ly2, maxR, maxS, mBS, minimz, minmax,
     &     MnrHdP, MnrHdS, MnrPrt, nDegen, nFac, ngQP0, nkx,
     &     nnb, nnCon0, nnCon, nnObj, nnJac, numLC, numLIQ,
     &     PreCon, Prob, QPmode, QPslvr, qpStat, sqStat
      double precision
     &     degen, dnormi, ObjLP, tolFP, tolQP, tolx, wtInf0,
     &     piNorm, pNorm1, pNorm2, rgNorm, sclObj, vimax, xNorm, Fx(1)
      external
     &     dnormi
!     ------------------------------------------------------------------
      integer            QPChol,     CG,         QN
      parameter         (QPChol = 0, CG     = 1, QN   = 2)
      integer            FP,         LP,         QP
      parameter         (FP     = 0, LP     = 1, QP   = 2)
      integer            Fix
      parameter         (Fix    = 0)
      integer            NO,         YES
      parameter         (NO     = 0, YES    = 1)
      integer            SaveB,      PrintS,     Wrap
      parameter         (SaveB  = 0, PrintS = 1, Wrap = 1)
      integer            Stats
      parameter        ( Stats  = 1 )

      parameter         (lvlScl =  75) ! scale option
      parameter         (eigH   = 200) ! type of QP Hessian
      parameter         (QPmode = 208) ! Current QP solver
      parameter         (PreCon = 209) ! Current precon mode
      parameter         (nFac   = 210) ! # of LU factorizations
      parameter         (linesL = 220) ! # lines in log     file
      parameter         (linesS = 221) ! # lines in summary file
      parameter         (MnrHdP = 223) ! >0  => Mnr heading for iPrint
      parameter         (MnrHdS = 225) ! >0 => Minor heading for iSumm
      parameter         (gotFac = 230) ! Save the LU factors
      parameter         (gotHes = 231) ! Save the reduced Hessian
      parameter         (qpStat = 235) ! QP user-routine call-status
      parameter         (cgItns = 386) ! Number of symmlq iterations
      parameter         (cgItn  = 387) ! symmlq itns for last minor

      double precision   zero,            one
      parameter         (zero   = 0.0d+0, one = 1.0d+0)
!     ------------------------------------------------------------------
      tolFP     = rw( 51) ! Minor Phase 1 Opt tol
      tolQP     = rw( 52) ! Minor Phase 2 Opt tol
      tolx      = rw( 56) ! Minor feasibility tolerance.
      wtInf0    = rw( 88) ! infeasibility weight

      nnObj     = iw( 22) ! # of objective variables
      lenR      = iw( 28) ! R(lenR) is the reduced Hessian factor
      maxR      = iw( 52) ! max columns of R.
      maxS      = iw( 53) ! max # of superbasics
      QPslvr    = iw( 55) ! = 0:1:2   => QPChol:CG:QN QP solver
      lEmode    = iw( 56) ! >0    => use elastic mode
      lvlInf    = iw( 73) ! Elastic option
      minmax    = iw( 87) ! 1, 0, -1  => MIN, FP, MAX
      itnlim    = iw( 89) ! limit on total iterations
      MnrPrt    = iw( 93) ! Minor print level
      iNewB     = iw(124) ! new basis file
      minimz    = iw(199) ! 1 (-1)    => minimize (maximize)
      nkx       = iw(247) ! dimension of kx and its inverse, kxN

      mProb     = cw( 51) ! Problem name

      ! Addresses

      lkx       = iw(251) ! j  = kx (jN) => col j of Jcol is variable jN
      lhfeas    = iw(284) ! hfeas(mBS)  = feasibility types
      lhEsta    = iw(285) ! hEstat(nb)  = status of elastics
      lkBS      = iw(292) ! kBS(mBS)    = ( B  S ) list
      lblBS     = iw(273) ! blBS(mBS)   = lower bounds for xBS
      lbuBS     = iw(274) ! buBS(mBS)   = upper bounds for xBS
      lpBS      = iw(277) ! pBS(nb)     = search direction
      lxBS      = iw(301) ! xBS(mBS)    = basics, superbasics
      lgQP      = iw(290) ! gQP(ngQP)   = QP gradient
      lgBS      = iw(291) ! gBS(mBS)    = BS components of g
      lrg       = iw(293) ! rg (maxS)   =  reduced gradient
      lrg2      = iw(294) ! rg2(maxS)   =  reduced gradient copy
      lR        = iw(295) ! R(lenR)     = factor of Z'HZ
      lAscal    = iw(296) ! Ascale(nb)  = row and column scales
      liy       = iw(308) ! iy (nb)     =  integer work vector
      liy1      = iw(309) ! iy1(nb)     =  integer work vector
      liy2      = iw(310) ! iy2(nb)     =  integer work vector
      ly        = iw(311) !  y (nb)     =  real work vector
      ly1       = iw(312) !  y1(nb)     =  real work vector
      ly2       = iw(313) !  y2(nb)     =  real work vector
      lHdx      = iw(288) ! Hdx(nnH)    = product of H with  x - x0
      lblSav    = iw(275) ! blSav(m)    = temp bounds
      lbuSav    = iw(276) ! buSav(m)    = temp bounds

      iExit = 0
      mBS   = m + maxS

      ! Figure out what type of problem we have.

      if (minmax .eq. 0  .or.(lEmode .eq. 2  .and.  lvlInf .eq. 2)) then
         Prob = FP
      else if (ngQP .eq. 0) then ! No explicit objective. Must be an LP.
         if (iObj .eq. 0) then
            Prob = FP
         else
            Prob = LP
         end if
      else !  Explicit objective. Check for quadratic term.
         if (nnH .gt. 0) then
            Prob = QP
         else
            Prob = LP
         end if
      end if

      iw(MnrHdP) = 0            ! Print the header for the Print   file
      iw(MnrHdS) = 0            ! Print the header for the summary file
      iw(linesL) = 0            ! Line count for the print   file
      iw(linesS) = 0            ! Line count for the summary file

      ! Initialize counters based on gotHes and gotFac (set in s3prtQ)

      if (iw(gotFac) .le. 0) then
         iw(nFac) = 0
      end if

      if (iw(gotHes) .le. 0) then
         Hcalls   = 0
      end if

      iw(cgItns) = 0
      iw(cgItn ) = 0

      itn        = 0
      itQP       = 0
      itQPmax    = itnlim
      nDegen     = 0
      nnCon      = 0
      nnCon0     = 1
      nnJac      = 0
      numLC      = m
      ngQP0      = max( ngQP , 1 )

      iw(eigH)   = 0            ! QP Hessian may or may not be definite
      iw(QPmode) = QPslvr       ! Local value of QPslvr
      ObjQP      = zero
      sclObj     = one

      ! Initialize quantities to avoid them being used before being set.

      call dload
     &   ( m    , zero, pi        , 1 )
      call dload
     &   ( ngQP0, zero, rw(lgQP)  , 1 )
      call iload
     &   ( nb   ,    0, iw(lhEsta), 1 )

      !-----------------------------------------------------------------
      ! Print the matrix statistics.
      ! Find the rowtypes for use in s5getB (they are held in iy2).
      !-----------------------------------------------------------------
      call s2Amat
     &   ( Stats, MnrPrt, m, n, nb,
     &     nnCon, nnJac, nnObj, iObj, numLC, numLIQ,
     &     ne, nlocA, locA, indA, Acol,
     &     bl, bu, iw(liy2),
     &     iw, leniw, rw, lenrw )

      !=================================================================
      ! Find a basis kBS(1:m) for the linear constraints and bounds.
      !=================================================================
      ! s5getB does the following.
      !  1. The linear constraints are (optionally) scaled.
      !  2. Elements x(n+1:n+m) of the initial x are assigned.
      !  3. An LP is used to find a feasible x for the bounds and
      !     linear equality constraints.
      !  The base point x0 is not touched.

      call s5getB
     &   ( inform, Start, QPlog, needB, m, maxS, mBS,
     &     n, nb, nnCon, nnJac, nnObj, nName, nS, itQP, itQPmax, itn,
     &     nDegen, numLC, numLIQ, tolFP, tolQP, tolx,
     &     nInf, sInf, wtInf0, iObj, sclObj, piNorm, rgNorm,
     &     ne, nlocA, locA, indA, Acol,
     &     hEtype, iw(lhEsta), iw(liy2), iw(lhfeas), hs, iw(lkBS),Names,
     &     rw(lAscal), bl, bu,rw(lblBS),rw(lbuBS),rw(lblSav),rw(lbuSav),
     &     rw(lgBS), pi, rc, nrhs0, nrhs, rhs,
     &     1, 0, x0, x, rw(lxBS),
     &     iw(liy), iw(liy1), rw(ly), rw(ly1), rw(ly2),
     &     cw, lencw, iw, leniw, rw, lenrw )

      if (ngObj .gt. 0  .and.  iw(lvlScl) .gt. 0) then
         call ddscl
     &      ( ngObj, rw(lAscal), 1, gObj, 1 )
      end if

!     Possible inform values are = -3,-2,-1, 0, >0

      if (inform .ne. 0) then
         if (inform .gt. 0) then
            iExit = inform      ! fatal error
         else if (inform .eq. -3) then
            iExit = 31          ! too many iterations
         else
            iExit = 12          ! infeasible linear equalities
         end if
      end if

      if (iExit .ne. 0) go to 900

!     ==================================================================
!     Solve the problem.
!     ==================================================================
      call s1time
     &   ( 2, 0, iw, leniw, rw, lenrw )

      useQP = (ngObj .gt. 0 .and. Prob .eq. LP)  .or.  Prob .eq. QP

      if (useQP) then
         if (   iw(QPmode) .eq. QPChol) then
            call s5sQP
     &         ( inform, Hprod, Hprod1, QPlog, gotR, Prob,
     &           lenR, m, maxS, mBS, n, nb, nDegen, Hcalls,
     &           ngQP0, ngQP, ngObj0, ngObj, nnH0, nnH, nS,
     &           itQP, itQPmax, itn,
     &           minimz, iObj, sclObj, ObjAdd, ObjQP,
     &           tolFP, tolQP, tolx, nInf, sInf, wtInf0, piNorm,
     &           ne, nlocA, locA, indA, Acol,
     &           hEtype, iw(lhEsta), iw(lhfeas), hs, iw(lkBS),
     &           rw(lAscal), bl, bu, rw(lblBS), rw(lbuBS),
     &           rw(lgBS), gObj, rw(lgQP), rw(lHdx), rw(lpBS),
     &           pi, rw(lR), rc, rw(lrg),
     &           nrhs0, nrhs, rhs, lenx0, nx0, x0,
     &           x, rw(lxBS), x, ! xFreez = x, not used
     &           iw(liy), iw(liy1), rw(ly), rw(ly1), rw(ly2),
     &           cu, lencu, iu, leniu, ru, lenru,
     &           cw, lencw, iw, leniw, rw, lenrw )
         else if (iw(QPmode) .eq. QN) then
            gotR = iw(gotHes) .gt. 0
            call s5sQN
     &         ( inform, Hprod, Hprod1, QPlog, gotR, Prob,
     &           lenR, m, maxS, mBS, n, nb, nDegen, Hcalls,
     &           ngQP0, ngQP, ngObj0, ngObj, nnH0, nnH, nS,
     &           itQP, itQPmax, itn,
     &           minimz, iObj, sclObj, ObjAdd, ObjQP,
     &           tolFP, tolQP, tolx, nInf, sInf, wtInf0, piNorm,
     &           ne , nlocA, locA, indA, Acol,
     &           hEtype, iw(lhEsta), iw(lhfeas), hs, iw(lkBS),
     &           rw(lAscal), bl, bu, rw(lblBS), rw(lbuBS),
     &           rw(lgBS), gObj, rw(lgQP), rw(lHdx), rw(lpBS),
     &           pi, rw(lR), rc, rw(lrg), rw(lrg2),
     &           nrhs0, nrhs, rhs, lenx0, nx0, x0,
     &           x, rw(lxBS), x, ! xFreez = x, not used
     &           iw(liy), iw(liy1), rw(ly), rw(ly1), rw(ly2),
     &           cu, lencu, iu, leniu, ru, lenru,
     &           cw, lencw, iw, leniw, rw, lenrw )
         end if

         if (inform .eq. 33  .and.  maxR .lt. maxS) then
            if (     iw(QPmode) .eq. QPChol) then
               iw(PreCon) = NO  ! with no preconditioning
               gotR       = .false.
            else if (iw(QPmode) .eq. QN    ) then
               iw(PreCon) = YES ! with QN preconditioning
            end if
            iw(QPmode)    = CG  ! Switch to CG
         end if

         if (iw(QPmode) .eq. CG) then
            if (QPslvr  .eq. CG) gotR = .false.
            call s5sQN
     &         ( inform, Hprod, Hprod1, QPlog, gotR, Prob,
     &           lenR, m, maxS, mBS, n, nb, nDegen, Hcalls,
     &           ngQP0, ngQP, ngObj0, ngObj, nnH0, nnH, nS,
     &           itQP, itQPmax, itn,
     &           minimz, iObj, sclObj, ObjAdd, ObjQP,
     &           tolFP, tolQP, tolx, nInf, sInf, wtInf0, piNorm,
     &           ne , nlocA, locA, indA, Acol,
     &           hEtype, iw(lhEsta), iw(lhfeas), hs, iw(lkBS),
     &           rw(lAscal), bl, bu, rw(lblBS), rw(lbuBS),
     &           rw(lgBS), gObj, rw(lgQP), rw(lHdx), rw(lpBS),
     &           pi, rw(lR), rc, rw(lrg), rw(lrg2),
     &           nrhs0, nrhs, rhs, lenx0, nx0, x0,
     &           x, rw(lxBS), x, ! xFreez = x, not used
     &           iw(liy), iw(liy1), rw(ly), rw(ly1), rw(ly2),
     &           cu, lencu, iu, leniu, ru, lenru,
     &           cw, lencw, iw, leniw, rw, lenrw )
         end if
         iExit = inform

      else
         call s5fixS
     &      ( Fix, m, maxS, mBS, n, nb, nS, hs, iw(lkBS),
     &        bl, bu, rw(lblBS), rw(lbuBS), x, rw(lxBS) )
         call s5sLP
     &      ( iExit, QPlog, Prob,
     &        m, mBS, n, nb, nDegen, itQP, itQPmax, itn,
     &        minimz, iObj, sclObj, ObjAdd,
     &        tolFP, tolQP, tolx, nInf, sInf, wtInf0, piNorm,
     &        ne, nlocA, locA, indA, Acol,
     &        hEtype, iw(lhEsta), iw(lhfeas), hs, iw(lkBS),
     &        rw(lAscal), bl, bu, rw(lblBS), rw(lbuBS),
     &        rw(lgBS), pi, rc,
     &        nrhs0, nrhs, rhs,
     &        x, rw(lxBS), x,   ! xFreez = x, not used
     &        iw(liy), iw(liy1), rw(ly), rw(ly1),
     &        cw, lencw, iw, leniw, rw, lenrw )

      end if
      call s1time
     &   (-2, 0, iw, leniw, rw, lenrw )

!     ==================================================================
!     Exit.
!     Set output variables and print a summary of the final solution.
!     ObjTru is printed in s4newB
!     ==================================================================
  900 call snWRAP( iExit, Solver, str, str2, iw, leniw )

      degen  = 100.0d+0 * nDegen / max( itn, 1 )

      ObjLP  = zero
      if (Prob .eq. FP) then
         ObjTru = zero
      else if (Prob .eq. LP  .or.  Prob .eq. QP) then
         ObjTru = ObjAdd

         if (iObj .gt. 0) then
            ObjLP  = x(n+iObj)*sclObj
            ObjTru = ObjTru + ObjLP
         end if

         if (ngQP .gt. 0) then
            ObjTru = ObjTru + ObjQP
         end if
      end if

      infsbl = nInf .gt. 0
      xNorm  = dnormi( n, x, 1 )

      ! Count basic nonlinear variables (used only for printing).

      nnb  = 0
      do j = 1, nnH
         if (hs(j) .eq. 3) nnb = nnb + 1
      end do

      if (inewB .gt. 0  .and.  iExit/10 .lt. 8) then
         k      = 1 + iExit/10
         call s4stat
     &      ( k, istate )
         call s4newB
     &      ( Wrap, iNewB, minimz, m, n, nb,
     &        nS, mBS, itn, nInf, sInf, ObjTru, iw(lkBS), hs,
     &        rw(lAscal), bl, bu, x, rw(lxBS), istate,
     &        cw, lencw, iw, leniw )
      end if

      ! Print statistics.

      call snPRNT(13,
     &     ' Problem name                 '//mProb, iw, leniw )
      write(str, 1900) itn, ObjTru
      call snPRNT( 3, str, iw, leniw )
      if (infsbl) then
         write(str, 1910) nInf, sInf
         call snPRNT( 3, str, iw, leniw )
      end if
      if (Prob .eq. QP) then
         write(str, 1920) Hcalls, ObjLP
         call snPRNT( 3, str, iw, leniw )
         write(str, 1930) ObjQP
         call snPRNT( 3, str, iw, leniw )
      end if
      if (nS .gt. 0) then
         write(str, 1970) nS, nnb
         call snPRNT( 3, str, iw, leniw )
      end if
      write(str, 1975) nDegen, degen
      call snPRNT( 3, str, iw, leniw )

!     ------------------------------------------------------------------
!     Unscale, save basis files and prepare to print the solution.
!     Clock 3 is "Output time".
!     ------------------------------------------------------------------
      call s1time
     &   ( 3, 0, iw, leniw, rw, lenrw )
      call s4savB
     &   ( inform, SaveB, minimz, m, n, nb, nkx,
     &     nnCon0, nnCon, ngQP0, ngQP, nName, nS,
     &     itn, nInf, sInf, wtInf0, vimax, iObj, sclObj, ObjTru,
     &     pNorm1, pNorm2, piNorm, xNorm,
     &     ne, nlocA, locA, indA, Acol, iw(lkx),
     &     iw(lhEsta), hs, rw(lAscal), bl, bu, Fx, rw(lgQP),
     &     Names, pi, rc, x,
     &     cw, lencw, iw, leniw, rw, lenrw )

      if (ngObj .gt. 0  .and.  iw(lvlScl) .gt. 0)
     &     call dddiv ( ngObj, rw(lAscal), 1, gObj, 1 )

!     If task = PrintS, s5savB prints the solution under the control
!     of lprSol (set by the  Solution  keyword in the SPECS file).
!     The printed solution may or may not be wanted, as follows:
!
!     lprSol = 0   means      No
!            = 1   means      If optimal, infeasible or unbounded
!            = 2   means      Yes
!            = 3   means      If error condition

      call s4savB
     &   ( inform, PrintS, minimz, m, n, nb, nkx,
     &     nnCon0, nnCon, ngQP0, ngQP, nName, nS,
     &     itn, nInf, sInf, wtInf0, vimax, iObj, sclObj, ObjTru,
     &     pNorm1, pNorm2, piNorm, xNorm,
     &     ne, nlocA, locA, indA, Acol, iw(lkx),
     &     iw(lhEsta), hs, rw(lAscal), bl, bu, Fx, rw(lgQP),
     &     Names, pi, rc, x,
     &     cw, lencw, iw, leniw, rw, lenrw )
      call s1time
     &   (-3, 0, iw, leniw, rw, lenrw )

!     ------------------------------------------------------------------
!     Set Obj for output.
!     Call  Hx  one last time with  nState .ge. 2.
!     Everything has been  unscaled, so we have to disable scaling.
!     ------------------------------------------------------------------
      lsSave     = iw(lvlScl)
      iw(lvlScl) = 0
      sqStat     = 2 + min( iExit/10,4 )
      iw(qpStat) = sqStat

      ObjLP = 0
      if (Prob .eq. FP) then
         ObjTru = zero
      else if (Prob .eq. LP  .or.  Prob .eq. QP) then
         ObjTru = ObjAdd

         if (iObj .gt. 0) then
            ObjLP  = x(n+iObj)*sclObj
            ObjTru = ObjTru + ObjLP
         end if

         if (ngQP .gt. 0) then
            call s5QPfg
     &         ( Hprod, Hprod1,
     &           ngQP, ngObj0, ngObj, nnH, sqStat, Hcalls, ObjQP,
     &           gObj, rw(lgQP), lenx0, nx0, x0, x, rw(ly),
     &           cu, lencu, iu, leniu, ru, lenru,
     &           cw, lencw, iw, leniw, rw, lenrw )
            ObjTru = ObjTru + ObjQP
         end if
      end if

      iw(lvlScl) = lsSave

!     Save some things needed by solvers calling SQOPT

  999 rw(421) = ObjTru          ! The true objective
      rw(422) = piNorm          ! Lagrange multiplier norm
      rw(423) = xNorm

      iw(421) = itn             ! Total iteration count
      iw(423) = maxS            ! max # of superbasics

      return

 1900 format(
     &     ' No. of iterations', i20, 2x,
     &     ' Objective value', 1p, e22.10)
 1910 format(
     &     ' No. of infeasibilities', i15, 2x,
     &     ' Sum of infeas', 1p, e24.10)
 1920 format(
     &     ' No. of Hessian products', i14, 2x,
     &     ' Objective row', 3x, 1p, e21.10)
 1930 format(
     &       40x,
     &     ' Quadratic objective', 1p, e18.10)
 1970 format(
     &     ' No. of superbasics', i19, 2x,
     &     ' No. of basic nonlinears', i14)
 1975 format(
     &     ' No. of degenerate steps', i14, 2x,
     &     ' Percentage', f27.2)

      end ! subroutine s5solv

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s5dflt
     &   ( m, n, lencObj, ncolH, cw, lencw, iw, leniw, rw, lenrw )

      implicit
     &     none
      integer
     &     lencObj, lencw, leniw, lenrw, m, n, ncolH, iw(leniw)
      double precision
     &     rw(lenrw)
      character
     &     cw(lencw)*8

!     ==================================================================
!     s5dflt checks and possibly prints the optional parameter values
!     for sqopt.
!
!     Optional parameters are checked and, if necessary,  changed to
!     reasonable values.
!
!     Note that parameters are checked before the amount of working
!     storage has been defined.
!
!     See  snworkspace.info  for full documentation of cw, iw and rw.
!
!     15 Nov 1991: first version.
!     02 Aug 2003: snPRNT adopted.
!     22 Jun 2004: Added default LU mod singularity tol
!     21 Dec 2004: Default LU tols fixed up.
!     02 May 2006: lvlTim removed.
!     01 Sep 2007: stkyOp added.
!     18 Jan 2010: PreCon initialized.
!     ==================================================================
      logical
     &     linear, QP
      integer
     &     cgItmx, iCrash, iBack, iDump, iLoadB, iMPS, iNewB, iInsrt,
     &     iOldB, iPnch, iPrint, iReprt, iSoln, itnlim, kchk,
     &     kDegen, kFac, klog, kReset, ksav, kSumm, lEmode, lprDbg,
     &     lprPrm, lprScl, lprSol, LUprnt, lvlInf, lvlPre, lvlPiv,
     &     lvlScl, lvlSys, maxmn, maxCol, maxR, maxS, mflush,
     &     minimz, minmax, minPrc, mMinor, MjrPrt, MnrPrt, mSkip,
     &     mNewSB, never, nout, nParPr, nPr1,
     &     nPr2, Precon, QPslvr, stkyOp, TPivot
      double precision
     &     infBnd, bigdx, bigFx, c4, c6, chzbnd, Dens1, Dens2,
     &     Lmax1, Lmax2, eps, eps0, eps1, eps2, eps3, eps4, etarg,
     &     Hcndbd, rmaxS, scltol, small, tCrash, tolCG, tolCon,
     &     tolDcp, tolDdp, tolDpp, tolDrp, tolDup, toldj3,
     &     tolFac, tolFP, tolNLP, tolpiv, tolQP, tolRow, tolSwp, tolUpd,
     &     tolx, Uspace, Utol1, Utol1m, Utol2, Utol2m, wtInf0, xdlim,
     &     Zcndbd
!     ------------------------------------------------------------------
      integer            QPChol,     CG
      parameter         (QPChol = 0, CG = 1)
      integer            idummy
      parameter         (idummy = -11111)
      double precision   zero,             one
      parameter         (zero   =  0.0d+0, one    = 1.0d+0)
      double precision   ten
      parameter         (ten    = 10.0d+0)
      double precision   tenp6,            hundrd
      parameter         (tenp6  = 1.0d+6,  hundrd = 100.0d+0)
!     ------------------------------------------------------------------
!     Set some local machine-dependent constants.

      eps       = rw(  1) ! machine precision.  IEEE DP  2.22e-16
      eps0      = rw(  2) ! eps**(4/5)          IEEE DP  3.00e-13
      eps1      = rw(  3) ! eps**(2/3)          IEEE DP  3.67e-11
      eps2      = rw(  4) ! eps**(1/2)          IEEE DP  1.49e-08
      eps3      = rw(  5) ! eps**(1/3)          IEEE DP  6.05e-06
      eps4      = rw(  6) ! eps**(1/4)          IEEE DP  1.22e-04

!     ------------------------------------------------------------------
!     rw(51)--rw(150): optional parameters set via the specs file.
!     ------------------------------------------------------------------
      tolFP     = rw( 51) ! Minor Phase 1 Opt tol
      tolQP     = rw( 52) ! Minor Phase 2 Opt tol
      tolCG     = rw( 54) ! cg tolerance

      tolx      = rw( 56) ! Minor feasibility tolerance

      tolpiv    = rw( 60) ! excludes small elements of y
      tolrow    = rw( 61) ! tolerance for the row error
      tCrash    = rw( 62) ! crash tolerance
      Utol1m    = rw( 63) ! abs tol for small diag of U in LU mod
      Utol2m    = rw( 64) ! rel tol for small diag of U in LU mod
      tolswp    = rw( 65) ! LU swap tolerance
      tolFac    = rw( 66) ! LU factor tolerance
      tolUpd    = rw( 67) ! LU update tolerance
      infBnd    = rw( 70) ! definition of plus infinity
      bigFx     = rw( 71) ! unbounded objective
      bigdx     = rw( 72) ! unbounded step
      lvlPre    = iw( 77) ! >0    => QN preconditioned CG
      xdlim     = rw( 80) ! Step limit
      etarg     = rw( 83) ! Quasi-Newton QP rg tolerance
      Hcndbd    = rw( 85) ! bound on the condition of Hz
      Zcndbd    = rw( 86) ! bound on the condition of Z
      wtInf0    = rw( 88) ! infeasibility weight

      scltol    = rw( 92) ! scale tolerance.
!     ------------------------------------------------------------------
!     rw(151)--rw(180) contain  parmLU  parameters for LUSOL.
!     ------------------------------------------------------------------
      Lmax1     = rw(151) ! max L-multiplier in factor
      Lmax2     = rw(152) ! max L-multiplier in update
      small     = rw(153) ! defn of small real
      Utol1     = rw(154) ! abs tol for small diag of U
      Utol2     = rw(155) ! rel tol for small diag of U
      Uspace    = rw(156) ! limit on waste space in U
      Dens1     = rw(157) ! switch to search maxcol columns and no rows
      Dens2     = rw(158) ! switch to dense LU
!     ------------------------------------------------------------------
!     rw(181)--rw(199) pass parameters into various routines.
!     ------------------------------------------------------------------
!     toldj3    = rw(186) ! current optimality tol
!     ------------------------------------------------------------------
!     iw(1)--iw(50): I/O file numbers and dimensions.
!     ------------------------------------------------------------------
      iPrint    = iw( 12) ! Print file
!     ------------------------------------------------------------------
!     iw(51)--iw(150): optional parameters set via the specs file.
!     ------------------------------------------------------------------
      maxR      = iw( 52) ! max columns of R.
      maxS      = iw( 53) ! max # of superbasics
      QPslvr    = iw( 55) ! 0(1) => QP(QN) QP solver
      lEmode    = iw( 56) ! >0    => use elastic mode
      kchk      = iw( 58) ! check (row) frequency
      kFac      = iw( 59) ! factorization frequency
      ksav      = iw( 60) ! save basis map
      klog      = iw( 61) ! log/print frequency
      kSumm     = iw( 62) ! Summary print frequency
      kDegen    = iw( 63) ! max. expansions of featol
      kReset    = iw( 64) ! Hessian frequency
      mFlush    = iw( 66) ! Hessian flush
      mSkip     = iw( 67) ! # largest value of nSkip
!     lvlSrt    = iw( 69) ! = 0:1:2:3 => cold:warm:basis:hot start
      lvlSys    = iw( 71) ! > 0   => print system info
      lvlInf    = iw( 73) ! Elastic option
      lvlScl    = iw( 75) ! scale option
      lvlPiv    = iw( 80) ! 0(1) LU threshold partial(complete) pivoting
      lprPrm    = iw( 81) ! > 0    => parms are printed
      lprScl    = iw( 83) ! > 0    => print the scales
      lprSol    = iw( 84) ! > 0    => print the solution
      lprDbg    = iw( 85) ! > 0    => private debug print
      minmax    = iw( 87) ! 1, 0, -1  => MIN, FP, MAX
      iCrash    = iw( 88) ! Crash option
      itnlim    = iw( 89) ! limit on total iterations
      mMinor    = iw( 91) ! limit on minor iterations
      MnrPrt    = iw( 93) ! Minor print level
      nParPr    = iw( 94) ! # of partial pricing sections
      mNewSB    = iw( 95) ! # of working set changes
      cgItmx    = iw( 97) ! CG iteration limit
      stkyOp    = iw(116) ! > 0 => optional parameters are sticky
      iBack     = iw(120) ! backup file
      iDump     = iw(121) ! dump file
      iLoadB    = iw(122) ! load file
      iMPS      = iw(123) ! MPS file
      iNewB     = iw(124) ! new basis file
      iInsrt    = iw(125) ! insert file
      iOldB     = iw(126) ! old basis file
      iPnch     = iw(127) ! punch file
      iReprt    = iw(130) ! Report file
      iSoln     = iw(131) ! Solution file
!     ------------------------------------------------------------------
!     iw(151)--iw(180) contain luparm parameters for LUSOL.
!     ------------------------------------------------------------------
      nout      = iw(151) ! unit # for printed messages
      LUprnt    = iw(152) ! print level in LU routines
      maxcol    = iw(153) ! lu1fac: max. # columns
!     ------------------------------------------------------------------
      PreCon    = iw(209) ! Current precon mode (based on QPslvr)

      c4        = max( 1.0d-4, eps3 )
      c6        = max( 1.0d-6, eps2 )
      never     = 99999999
      QP        = ncolH .gt. 0
      linear    = .not. QP

!     ==================================================================
!     Check the optional parameters.
!     ==================================================================
      if (iBack  .eq. idummy ) iBack  =     0
      if (iDump  .eq. idummy ) iDump  =     0
      if (iLoadB .eq. idummy ) iLoadB =     0
      if (iNewB  .eq. idummy ) iNewB  =     0
      if (iInsrt .eq. idummy ) iInsrt =     0
      if (iOldB  .eq. idummy ) iOldB  =     0
      if (iPnch  .eq. idummy ) iPnch  =     0
      if (iReprt .eq. idummy ) iReprt =     0
      if (iSoln  .eq. idummy ) iSoln  =     0

!     Set unspecified frequencies or silly values to defaults.

      if (kchk   .eq. idummy ) kchk   =    60
      if (kfac   .le.    0   ) then
                               kfac   =   100
                     if ( QP ) kfac   =    50
      end if
      if (klog  .eq. idummy  ) klog   =   100
      if (kSumm .eq. idummy  ) kSumm  =   100
      if (ksav  .eq. idummy  ) ksav   =   100
      if (kDegen.eq. idummy  ) kDegen = 10000

!     Sometimes, frequency 0 means "almost never".

      if (kchk   .le. 0      ) kchk   = never
      if (klog   .le. 0      ) klog   = never
      if (ksav   .le. 0      ) ksav   = never
      if (kSumm  .le. 0      ) kSumm  = never
      if (kDegen .le. 0      ) kDegen = never
      if (iCrash .lt. 0      ) iCrash =  3
      if (minmax .eq. idummy ) minmax =  1
      if (minmax .eq. -1     ) then
                               minimz = -1
      else
                               minimz =  1
      end if

      if (mMinor .lt. 0      ) mMinor = max(1000, 5*max(n,m))
      if (mNewSB .le. 0      ) mNewSB = never
      if (lprDbg .lt. 0      ) lprDbg = 0
      if (lprPrm .lt. 0      ) lprPrm = 1
      if (lprScl .lt. 0      ) lprScl = 0
      if (lprSol .lt. 0      ) lprSol = 2
!     lvlSrt is checked in s3argA or s3argB
!     if (lvlSrt .lt. 0      ) lvlSrt = 0
      if (MnrPrt .lt. 0      ) MnrPrt = 1
                               MjrPrt = MnrPrt
      if (lvlInf .lt. 0  .or.  lvlInf .gt. 2
     &                       ) lvlInf = idummy
      if (lvlInf .eq. idummy ) lvlInf = 1
      if (lvlSys .lt. 0      ) lvlSys = 0
      if (lEmode .lt. 0  .or.  lEmode .gt. 2
     &                       ) lEmode = idummy
      if (lEmode .eq. idummy ) lEmode = 1

      if (stkyOp .lt. 0      ) stkyOp = 0

!     Check superbasics limit and reduced Hessian size.

      if ( QP ) then
         if (maxR .lt. 0     ) maxR   = min( 2000, ncolH+1 )
         if (maxS .lt. 0     ) maxS   =            ncolH+1
                               maxR   = max( min( maxR ,n ) , 0 )
                               maxS   = max( min( maxS ,n ) , 1 )
      else ! linear
         if (maxS   .le. 0   ) maxS   = 1
         if (maxR   .le. 0   ) maxR   = 1
      end if

      if (maxS   .lt. maxR   ) maxS   = maxR

      if (QPslvr .lt. 0      ) QPslvr = QPChol
      if (maxR   .eq. 0      ) QPslvr = CG
      if (lvlPre .lt. 0 .or.
     &    lvlPre .gt. 1      ) then
                               lvlPre = 0
                               PreCon = 0
      else
                               PreCon = 1
      end if
      if (cgItmx .lt. 0      ) cgItmx = 100
      if (etarg  .lt. zero  .or.
     &    etarg  .gt. one    ) etarg  = 0.5d+0

!     Check other options.

      if (lvlScl .lt. 0      ) lvlScl = 2
                               lvlScl = min( lvlScl, 2 )

      if (nParPr .le. 0      ) nParPr = 10
                               minPrc = 10
                               nPr1   = n / nParPr
                               nPr2   = m / nParPr
      if (max( nPr1, nPr2 ) .lt. minPrc) then
                               maxmn  = max( m, n )
                               nParPr = maxmn / min( maxmn, minPrc )
      end if

      rmaxS  = maxS
      cHzbnd = max ( one/(hundrd*eps*rmaxS), tenp6 )

      if (infBnd   .lt. zero ) infBnd = 1.0d+20
      if (bigFx    .le. zero ) bigFx  = 1.0d+15
      if (bigdx    .le. zero ) bigdx  = infBnd
      if (Hcndbd   .le. zero ) Hcndbd = cHzbnd
      if (xdlim    .le. zero ) xdlim  = 2.0d+0
      if (Zcndbd   .le. zero ) then
          if (QPslvr .eq. QPChol) then
                               Zcndbd = 1.0d+4
          else
                               Zcndbd = 1.0d+6
          end if
      end if

      if (tCrash   .lt. zero  .or.
     &    tCrash   .ge. one  ) tCrash = 0.1d+0

!     ------------------------------------
!     Set up the parameters for lu1fac.
!     ------------------------------------
      if (maxcol .lt.  0     ) maxcol =   5
      if (LUprnt .eq.  idummy) LUprnt =  -1
                               nout   =  iPrint
      if (lvlSys .eq.  0     ) nout   =  0
      if (MnrPrt .gt. 10     ) LUprnt =  0
      if (lprDbg .eq. 51     ) LUprnt =  1
      if (lprDbg .eq. 52     ) LUprnt =  2
      if (iPrint .lt.  0     ) LUprnt = -1
      if (lvlPiv .le.  0     ) lvlPiv =  0
      if (lvlPiv .gt.  3     ) lvlPiv =  0
                               TPivot =  lvlPiv
      if (linear) then
                               tolDpp =  hundrd
                               tolDrp =  ten
                               tolDcp =  ten
                               tolDdp =  ten
                               tolDup =  ten
      else ! QP
                               tolDpp =  3.99d+0
                               tolDrp =  3.99d+0
                               tolDcp =  3.99d+0
                               tolDdp =  3.99d+0
                               tolDup =  3.99d+0
      end if
      if (tolFac .lt. one    ) then
         if (lvlPiv .eq.   0 ) tolFac =  tolDpp
         if (lvlPiv .eq.   1 ) tolFac =  tolDrp
         if (lvlPiv .eq.   2 ) tolFac =  tolDcp
         if (lvlPiv .eq.   3 ) tolFac =  tolDdp
      end if
      if (tolUpd    .lt. one ) tolUpd =  tolDup
                               Lmax1  =  tolFac
                               Lmax2  =  tolUpd
      if (Utol1     .le. zero) Utol1  =  eps1
      if (Utol2     .le. zero) Utol2  =  eps1
      if (Utol1m    .le. zero) Utol1m =  eps1
      if (Utol2m    .le. zero) Utol2m =  eps1
      if (Dens2     .lt. zero) Dens2  =  0.6d+0
      if (small     .le. zero) small  =  eps0
      if (Uspace    .le. zero) Uspace =  3.0d+0
      if (Dens1     .le. zero) Dens1  =  0.3d+0

!     Set some tolerances.
!     Set the optimality tolerance.
!     Solve the QP subproblems fairly accurately.

      if (tolCG  .le. zero   ) tolCG  =  1.0d-2
      if (tolQP  .le. zero   ) tolQP  =  c6
      if (tolFP  .le. zero   ) tolFP  =  tolQP
      if (tolrow .le. zero   ) tolrow =  c4
      if (tolswp .le. zero   ) tolswp =  eps4
      if (tolx   .le. zero   ) tolx   =  c6
                               toldj3 =  tolQP
      if (scltol .le. zero   ) scltol =  0.90d+0
      if (scltol .ge. one    ) scltol =  0.99d+0
      if (tolpiv .le. zero   ) tolpiv =  eps1

      if (wtInf0 .lt. zero   ) wtInf0 =  1.0d+0

      if (iBack  .eq. iNewB  ) iBack  = 0
      if (itnlim .lt. 0      ) itnlim = max(10000, 10*max(n,m))

!     Load tolerances used to mark variables during printing in s4SavB.

      tolNLP  = tolQP
      tolCon  = tolx

!     ------------------------------------------------------------------
!     Re-assign the options to their respective work arrays.
!     ------------------------------------------------------------------
      rw( 51) = tolFP
      rw( 52) = tolQP
      rw( 53) = tolNLP
      rw( 54) = tolCG
      rw( 56) = tolx
      rw( 57) = tolCon
      rw( 60) = tolpiv
      rw( 61) = tolrow
      rw( 62) = tCrash
      rw( 65) = tolswp
      rw( 66) = tolFac
      rw( 67) = tolUpd
      rw( 70) = infBnd
      rw( 71) = bigFx
      rw( 72) = bigdx
      rw( 80) = xdlim
      rw( 83) = etarg
      rw( 85) = Hcndbd
      rw( 86) = Zcndbd
      rw( 88) = wtInf0
      rw( 92) = scltol

      rw(151) = Lmax1  ! max L-multiplier in factor
      rw(152) = Lmax2  ! max L-multiplier in update
      rw(153) = small  ! defn of small real
      rw(154) = Utol1  ! abs tol for small diag of U
      rw(155) = Utol2  ! rel tol for small diag of U
      rw(156) = Uspace ! limit on waste space in U
      rw(157) = Dens1  ! switch to search maxcol columns and no rows
      rw(158) = Dens2  ! switch to dense LU
      rw(181) = tolDpp
      rw(182) = tolDcp
      rw(183) = tolDup
      rw(186) = toldj3
      rw(187) = tolDrp

      iw( 52) = maxR
      iw( 53) = maxS
      iw( 55) = QPslvr
      iw( 56) = lEmode
      iw( 58) = kchk
      iw( 59) = kFac
      iw( 60) = ksav
      iw( 61) = klog
      iw( 62) = kSumm
      iw( 63) = kDegen
      iw( 64) = kReset
      iw( 66) = mFlush
      iw( 67) = mSkip
!     iw( 69) = lvlSrt
      iw( 71) = lvlSys
      iw( 73) = lvlInf
      iw( 75) = lvlScl
      iw( 77) = lvlPre
      iw( 80) = lvlPiv
      iw( 81) = lprPrm
      iw( 83) = lprScl
      iw( 84) = lprSol
      iw( 85) = lprDbg
      iw( 87) = minmax
      iw( 88) = iCrash
      iw( 89) = itnlim
      iw( 91) = mMinor
      iw( 92) = MjrPrt
      iw( 93) = MnrPrt
      iw( 94) = nParPr
      iw( 95) = mNewSB
      iw( 97) = cgItmx
      iw(116) = stkyOp
      iw(120) = iBack
      iw(121) = iDump
      iw(122) = iLoadB
      iw(123) = iMPS
      iw(124) = iNewB
      iw(125) = iInsrt
      iw(126) = iOldB
      iw(127) = iPnch
      iw(130) = iReprt
      iw(131) = iSoln
      iw(151) = nout
      iw(152) = LUprnt
      iw(153) = maxcol
      iw(156) = TPivot
      iw(199) = minimz

      iw(209) = PreCon   ! not optional parameters, but set here.
      rw(186) = toldj3

      end ! subroutine s5dflt

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s5Map
     &   ( m, n, nkx, ngObj, nnH,
     &     lenR, maxR, maxS,
     &     nextcw, nextiw, nextrw, iw, leniw )

      implicit
     &     none
      integer
     &     m, n, nkx, ngObj, nnH, lenR, maxR, maxS,
     &     nextcw, nextiw, nextrw, leniw, iw(leniw)

!     ==================================================================
!     s5Map   allocates all array storage for sqopt,
!     using the values:
!        m    , n    , ne
!        maxS                    Set in s5dflt.
!        ngObj, nnH              Set from the argument list.
!        lenR                    Set in the calling program.
!
!     15 Nov 1991: First version based on Minos 5.4 routine m2core.
!     12 Nov 1994: Converted to integer and real storage.
!     06 Aug 1996: First min sum version.
!     14 Jul 1997: Thread-safe version.
!     01 May 1998: First version called by sqMem. This simplified
!                  version may slightly overestimate needed memory.
!     02 Aug 2003: snPRNT adopted.
!     13 May 2005: Bug fix: ly3 assigned to iw correctly
!     18 Jun 2008: Added space for iy2, pBS and rg2.
!     ==================================================================
      integer
     &     mBS, nb, lAscal, lblBS, lbuBS, lblSav, lbuSav, ldx, lgBS,
     &     lgQP, lHdx, lhEsta, lhfeas, liy, liy1, liy2, lkBS, lkx,
     &     lQPrhs, lpBS, lR, lr1, lr2, lrg, lrg2, ls1, ls2, ls3,
     &     lxBS, lxScal, ly, ly1, ly2, ly3, ngQP
!     ------------------------------------------------------------------
      ngQP    = max( ngObj, nnH )
      mBS     = m     + maxS
      nb      = n     + m

!     sqopt can use all of cw, iw and rw
!     except the first user workspace partitions.

      lkx    = nextiw
      lhfeas = lkx    + nkx
      lkBS   = lhfeas + mBS
      lhEsta = lkBS   + mBS
      liy    = lhEsta + nb
      liy1   = liy    + nb
      liy2   = liy1   + nb
      nextiw = liy2   + nb

!     Addresses for the double precision arrays.

      lAscal = nextrw
      ly     = lAscal + nb
      ly1    = ly     + nb
      ly2    = ly1    + nb
      if (maxR .lt. maxS) then  ! Define SYMMLQ workspace
         ly3   = ly2  + nb
         ls1   = ly3  + nb
         ls2   = ls1  + maxS
         ls3   = ls2  + maxS
         lr1   = ls3  + maxS
         lr2   = lr1  + maxS
         lblBS = lr2  + maxS
      else
         ly3   = ly2  + nb
         ls1   = ly3
         ls2   = ls1
         ls3   = ls2
         lr1   = ls3
         lr2   = lr1
         lblBS = lr2
      end if
      lbuBS  = lblBS  + mBS
      lxBS   = lbuBS  + mBS
      lxScal = lxBS   + mBS
      lHdx   = lxScal + nnH
      lpBS   = lHdx   + nnH
      lgQP   = lpBS   + nb
      lgBS   = lgQP   + ngQP
      lR     = lgBS   + mBS
      lrg    = lR     + lenR
      lrg2   = lrg    + maxS
      lblSav = lrg2   + maxS
      lbuSav = lblSav + nb
      lQPrhs = lbuSav + nb
      ldx    = lQPrhs + m
      nextrw = ldx    + ngQP

!     ---------------------------
!     Store the addresses in iw.
!     ---------------------------
      iw(251) = lkx

      iw(273) = lblBS
      iw(274) = lbuBS
      iw(275) = lblSav
      iw(276) = lbuSav
      iw(277) = lpBS
      iw(278) = lQPrhs

      iw(284) = lhfeas
      iw(285) = lhEsta

      iw(287) = ldx
      iw(288) = lHdx
      iw(290) = lgQP
      iw(291) = lgBS
      iw(292) = lkBS
      iw(293) = lrg
      iw(294) = lrg2
      iw(295) = lR
      iw(296) = lAscal

      iw(301) = lxBS
      iw(302) = lxScal

      iw(308) = liy
      iw(309) = liy1
      iw(310) = liy2
      iw(311) = ly
      iw(312) = ly1
      iw(313) = ly2
      iw(314) = ly3

      iw(353) = lr1
      iw(354) = lr2
      iw(355) = ls1
      iw(356) = ls2
      iw(357) = ls3

      end ! subroutine s5Map

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s5sLP
     &   ( iExit, QPlog, Prob,
     &     m, mBS, n, nb, nDegen, itQP, itQPmax, itn,
     &     minimz, iObj, sclObj, ObjAdd,
     &     tolFP, tolQP, tolx, nInf, sInf, wtInf, piNorm,
     &     ne, nlocA, locA, indA, Acol,
     &     hEtype, hEstat, hfeas, hs, kBS,
     &     Ascale, bl, bu, blBS, buBS,
     &     gBS, pi, rc,
     &     nrhs0, nrhs, rhs,
     &     x, xBS, xFreez,
     &     iy, iy1, y, y1,
     &     cw, lencw, iw, leniw, rw, lenrw )

      implicit
     &     none
      external
     &     QPlog
      integer
     &     iExit, iObj, itn, itQP, itQPmax, lencw, leniw, lenrw, m,
     &     mBS, minimz, n, nb, ne, nInf, nDegen, nlocA, nrhs0, nrhs,
     &     Prob, hfeas(mBS), kBS(mBS), hEtype(nb), hEstat(nb), hs(nb),
     &     locA(nlocA), indA(ne), iy(nb), iy1(nb), iw(leniw)
      double precision
     &     sclObj, ObjAdd, sInf, wtInf, piNorm, tolFP, tolQP, tolx,
     &     Acol(ne), Ascale(nb), bl(nb), bu(nb), blBS(mBS), buBS(mBS),
     &     gBS(mBS), pi(m), rc(nb), rhs(nrhs0), x(nb), xBS(mBS),
     &     xFreez(nb), y(nb), y1(nb), rw(lenrw)
      character
     &     cw(lencw)*8

      !=================================================================
      ! s5sLP   solves the current problem.  An initial basis is assumed
      ! to be specified by nS, hs, x and the superbasic parts of kBS.
      !
      ! 25 Oct 2003: First version of s5sLP based on s5SQP.
      ! 08 Mar 2004: gotFac implemented.
      ! 18 Jun 2008: xFreez added as argument.
      ! 29 Apr 2011: Some minor cosmetic changes.
      !=================================================================
      character
     &     ProbTag*20, str*120
      logical
     &     badLU, done, Elastc, found, FPonly, getOpt, infeas, itnBig,
     &     LUok, needLU, needx, newB, newLU, unbnded
      integer
     &     gotFac, lEmod0, lEmode, inform, LUreq, lvlIn0, lvlInf,
     &     MnrHdP, MnrHdS, MnrPrt, nnH, nS, nSwap,
     &     PrtLvl, subopt, typeLU
      double precision
     &     rgNorm
      !-----------------------------------------------------------------
      integer            BT
      parameter         (BT     = 3)
      integer            FP,         LP,         QP
      parameter         (FP     = 0, LP     = 1, QP = 2)
      integer            No
      parameter         (No     =-1)
      parameter         (MnrHdP = 223) ! >0 => Mnr heading for iPrint
      parameter         (MnrHdS = 225) ! >0 => Mnr heading for iSumm
      !-----------------------------------------------------------------
      lEmode    = iw( 56) ! >0    => use elastic mode
      lvlInf    = iw( 73) ! Elastic option
      MnrPrt    = iw( 93) ! Minor print level
      gotFac    = iw(230) ! >0 => Save the LU factors
      PrtLvl    = MnrPrt

      ProbTag   = 'linear constraints'

      call iload ( nb, 0, hEstat, 1 )

      Elastc    = lEmode .eq. 2
      getOpt    = Prob   .eq. LP   .or.  Prob .eq. QP
      FPonly    = Prob   .eq. FP

      nS        = 0             ! Local value
      nnH       = 0             ! Local value

      !-----------------------------------------------------------------
      ! Call s5LP with argument "FP" to find a feasible point.
      !-----------------------------------------------------------------
      if (Elastc) then
         ! Phase 2 will start in elastic mode with elastic objective
         ! determined by lvlInf. For FP, make the nonelastics feasible.

         lEmod0 = lEmode
         lvlIn0 = lvlInf
      else
         ! FP will make all constraints feasible.
         ! Enter elastic mode if infeasible and minimize the elastic
         ! infeasibilities.

         lEmod0 = 1
         lvlIn0 = 2
      end if

      subopt = No

      if (gotFac .eq. 0) then
         LUreq  = 1
      else
         LUreq  = 0
      end if
      typeLU = BT
      LUok   = .true.
      done   = .false.

*     ==================================================================
*+    while (.not. done  .and.  LUok) do
  500 if    (.not. done  .and.  LUok) then

         needLU = LUreq .gt. 0

         call s2Bfac
     &      ( iExit, typeLU, needLU, newLU, newB,
     &        iObj, itn, PrtLvl, LUreq,
     &        m, mBS, n, nb, nnH, nS, nSwap,
     &        ne, nlocA, locA, indA, Acol,
     &        kBS, hs, bl, bu, blBS, buBS,
     &        nrhs0, nrhs, rhs, x, xBS,
     &        iy, iy1, y, y1, iw, leniw, rw, lenrw )
         if (iExit .ne. 0) go to 900

         needx  = needLU

         call s5LP
     &      ( inform, FP, ProbTag, Elastc, subopt,
     &        QPlog, needLU, needx,
     &        m, n, nb, nDegen, itQP, itQPmax, itn,
     &        lEmod0, lvlIn0, PrtLvl,
     &        minimz, iObj, sclObj, ObjAdd, tolFP, tolQP, tolx,
     &        nInf, sInf, wtInf, piNorm, rgNorm,
     &        ne, nlocA, locA, indA, Acol,
     &        hEtype, hEstat, hfeas, hs, kBS,
     &        Ascale, bl, bu, blBS, buBS,
     &        gBS, pi, rc, nrhs0, nrhs, rhs, x, xBS, xFreez,
     &        iy, iy1, y, y1,
     &        cw, lencw, iw, leniw, rw, lenrw )

         ! Check for trouble in s5LP.  Here are the possibilities:

         badLU   = inform .gt.  0 ! Fatal LU error
         found   = inform .eq.  0 ! LP solution found
         infeas  = inform .eq. -1 ! LP is infeasible
         unbnded = inform .eq. -2 ! LP is unbounded
         itnBig  = inform .eq. -3 ! Too many iterations

         if (badLU) then
            iExit = inform      ! Fatal LU error
            go to 900
         end if

         done = found .or. infeas .or. itnBig

         if (.not. done) then
            !===========================================================
            ! Trouble.
            ! The phase 1 was unbounded, which can only occur if a bad
            ! basis has given a large search direction.
            !===========================================================
            ! Assume infeasible constraints.  Factor with tighter tols.

            infeas = .true.

            write(str, 1000) itn
            call snPRNT( 3, str, iw, leniw )
            call s2tryLU
     &         ( itn, 22, nS, LUreq, LUok, typeLU,
     &           iw, leniw, rw, lenrw )
         end if
         go to 500
      end if
!+    end while
!     ==================================================================

      if (.not. found                        ) go to 800 ! Itns or inf
      if (nInf   .gt. 0  .and.  lEmode .ne. 1) go to 800 ! Infeas

      if (getOpt) then
!        call dcopy ( nb, x, 1, xFreez, 1 ) ! Save feasible nonelastics

!        if (PrtLvl .ge. 1  .and.  PrtLvl .lt. 10) then
         if (PrtLvl .ge. 10) then
            if ( Elastc ) then
               write(str, 1200) itn
               call snPRNT( 3, str, iw, leniw )
            else
               write(str, 1300) itn
               call snPRNT( 3, str, iw, leniw )
            end if
         end if
      end if

      !=================================================================
      ! Get optimal.
      !=================================================================

      if (Elastc) then
         ! Relax
      else if (lEmode .eq. 1) then ! Elastc mode never needed
         lvlInf = 1
      end if

      LUreq  = 0
      typeLU = BT
      LUok   = .true.
      done   = FPonly

!     ==================================================================
!+    while (.not. done  .and.  LUok) do
  600 if    (.not. done  .and.  LUok) then

         newLU  = LUreq .gt. 0
         needx  = needLU

         !--------------------------------------------------------------
         ! LP with objective row in A.
         !--------------------------------------------------------------
         iw(MnrHdP) = 0
         iw(MnrHdS) = 0
         call s5LP
     &      ( inform, Prob, ProbTag, Elastc, subopt,
     &        QPlog, needLU, needx,
     &        m, n, nb, nDegen, itQP, itQPmax, itn,
     &        lEmode, lvlInf, PrtLvl,
     &        minimz, iObj, sclObj, ObjAdd, tolFP, tolQP, tolx,
     &        nInf, sInf, wtInf, piNorm, rgNorm,
     &        ne, nlocA, locA, indA, Acol,
     &        hEtype, hEstat, hfeas, hs, kBS,
     &        Ascale, bl, bu, blBS, buBS,
     &        gBS, pi, rc, nrhs0, nrhs, rhs, x, xBS, xFreez,
     &        iy, iy1, y, y1,
     &        cw, lencw, iw, leniw, rw, lenrw )

         ! Check for trouble in s5LP.  Here are the possibilities:

         badLU   = inform .gt.  0 ! Fatal LU error
         found   = inform .eq.  0 ! LP solution found
         infeas  = inform .eq. -1 ! LP is infeasible
         unbnded = inform .eq. -2 ! LP is unbounded
         itnBig  = inform .eq. -3 ! Too many iterations

         if (badLU) then
            iExit = inform
            go to 900
         end if

         done    = found  .or. unbnded .or. itnBig

         if (done) then
            !===========================================================
            ! Relax, we are finished.
            !===========================================================
         else
            !-----------------------------------------------------------
            ! The non-elastics are infeasible. This should not happen.
            ! Phase 1 has already found a feasible point for the
            ! nonelastics, so the basis must be ill-conditioned.
            ! Refactorize with tighter tols and restart at the known
            ! feasible point.  Reduce the feasibility tol to try and
            ! prevent repeats.
            !-----------------------------------------------------------
            write(str, 1100) itn
            call snPRNT( 3, str, iw, leniw )
            call s2tryLU
     &         ( itn, 22, nS, LUreq, LUok, typeLU,
     &           iw, leniw, rw, lenrw )
         end if

         go to 600
      end if
!+    end while
!     ==================================================================

  800 if (found) then
         if (nInf .eq. 0) then
            if (FPonly) then
               iExit =  2       ! Feasible
            else
               iExit =  1       ! optimal
            end if
         else if (lEmode .eq. 0) then
            iExit = 11          ! infeasible linear constraints
         else
            iExit = 14          ! infeasibilites minimized
         end if
      else if (infeas ) then
            iExit = 11          ! infeasible nonelastics
      else if (unbnded) then
            iExit = 21          ! unbounded
      else if (itnBig ) then
            iExit = 31          ! too many iterations
      end if

  900 return

 1000 format(' Itn', i7, ': Infeasible nonelastics in feasibility',
     &                   ' phase')
 1100 format(' Itn', i7, ': Infeasible nonelastics in LP optimality',
     &                   ' phase')
 1200 format(' Itn', i7, ': Feasible non-elastics')
 1300 format(' Itn', i7, ': Feasible constraints')

      end ! subroutine s5sLP

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s5sQP
     &   ( iExit, Hprod, Hprod1, QPlog, gotR, Prob,
     &     lenR, m, maxS, mBS, n, nb, nDegen, Hcalls,
     &     ngQP0, ngQP, ngObj0, ngObj, nnH0, nnH, nS,
     &     itQP, itQPmax, itn,
     &     minimz, iObj, sclObj, ObjAdd, ObjQP,
     &     tolFP, tolQP, tolx, nInf, sInf, wtInf, piNorm,
     &     ne , nlocA, locA, indA, Acol,
     &     hEtype, hEstat, hfeas, hs, kBS,
     &     Ascale, bl, bu, blBS, buBS,
     &     gBS, gObj, gQP, Hdx, pBS, pi, R, rc, rg,
     &     nrhs0, nrhs, rhs, lenx0, nx0, x0,
     &     x, xBS, xFreez,
     &     iy, iy1, y, y1, y2,
     &     cu, lencu, iu, leniu, ru, lenru,
     &     cw, lencw, iw, leniw, rw, lenrw )

      implicit
     &     none
      external
     &     Hprod, Hprod1, QPlog
      logical
     &     gotR
      integer
     &     Hcalls, iExit, iObj, itn, itQP, itQPmax, lenR, lenx0, lencu,
     &     leniu, lenru, lencw, leniw, lenrw, m, maxS, mBS, minimz,
     &     n, nb, ne, nInf, ngQP0, ngQP, nnH0, nnH, nS, nDegen,
     &     ngObj0, ngObj, nlocA, nrhs0, nrhs, nx0, Prob,
     &     hfeas(mBS), kBS(mBS), hEtype(nb), hEstat(nb), hs(nb),
     &     locA(nlocA), indA(ne), iy(nb), iy1(nb), iu(leniu), iw(leniw)
      double precision
     &     sclObj, ObjAdd, ObjQP, sInf, wtInf, piNorm,
     &     tolFP, tolQP, tolx,
     &     Acol(ne), Ascale(nb), bl(nb), bu(nb), blBS(mBS), buBS(mBS),
     &     gBS(mBS), gObj(ngObj0), gQP(ngQP0), Hdx(nnH0), pBS(nb),
     &     pi(m), R(lenR), rc(nb), rg(maxS), rhs(nrhs0),
     &     x(nb), x0(lenx0), xBS(mBS), xFreez(nb),
     &     y(nb), y1(nb), y2(nb), ru(lenru), rw(lenrw)
      character
     &     cu(lencu)*8, cw(lencw)*8

      !=================================================================
      ! s5sQP   solves the current problem.  An initial basis is assumed
      ! to be specified by nS, hs, x and the superbasic parts of kBS.
      ! In particular, there must be nS values hs(j) = 2, and the
      ! corresponding j's must be listed in kBS(m+1) thru kBS(m+nS).
      ! The ordering in kBS matches the reduced Hessian R (if any).
      !
      ! 05 Oct 1994: First version of s5sQP.
      ! 06 Aug 1996: Min Sum option added.
      ! 14 Jul 1997: Thread-safe version.
      ! 28 Apr 2002: Updated to reflect LU changes.
      ! 02 Aug 2003: snEXIT and snPRNT adopted.
      ! 08 Mar 2004: gotFac implemented.
      ! 16 May 2006: Explicit target itQP added.
      ! 27 Apr 2011: Finds indefinite Hessian with one BS factorize.
      ! 28 Apr 2011: Do one refactorize after detecting indefiniteness.
      !=================================================================
      character
     &     ProbTag*20, str*120
      logical
     &     badLU, badZg, badZ, done, Elastc, found, FPonly, getOpt,
     &     indefH, indefQP, infeas, itnBig, LUok, needLU, needx, newB,
     &     newLU, nSBig, unbnded, weakMin
      integer
     &     gotFac, inform, itQPtargt, lEmod0, lEmode, LUreq,
     &     lvlIn0, lvlInf, MnrHdP, MnrHdS, MnrPrt, nSwap, PrtLvl,
     &     subopt, typeLU, zngQP, znnH
      double precision
     &     eps, flmax, ObjFP, plInfy, rgNorm, targtH, targtZ,
     &     Hcndbd, Zcndbd
      !-----------------------------------------------------------------
      integer            BS        , BT
      parameter         (BS     = 2, BT     = 3)
      integer            FP,         LP,         QP
      parameter         (FP     = 0, LP     = 1, QP = 2)
      integer            No
      parameter         (No     =-1)

      parameter         (MnrHdP = 223) ! >0 => Mnr heading for iPrint
      parameter         (MnrHdS = 224) ! >0 => Mnr heading for iSumm

      double precision   zero,            one
      parameter         (zero   = 0.0d+0, one    = 1.0d+0)
      !-----------------------------------------------------------------
      lEmode    = iw( 56) ! >0    => use elastic mode
      lvlInf    = iw( 73) ! Elastic option
      MnrPrt    = iw( 93) ! Minor print level
      gotFac    = iw(230) ! >0 => Save the LU factors

      eps       = rw(  1) ! machine precision.  IEEE DP  2.22e-16
      flmax     = rw(  8) ! est. of the largest pos. real
      Hcndbd    = rw( 85) ! bound on the condition of Hz
      Zcndbd    = rw( 86) ! bound on the condition of Z

      PrtLvl    = MnrPrt
      subopt    = No

      plInfy    = flmax
      ObjQP     = zero
      targtH    = Hcndbd
      targtZ    = Zcndbd
      itQPtargt = itQPmax

      ProbTag   = 'linear constraints'

      getOpt    = Prob   .eq. LP   .or.  Prob .eq. QP
      FPonly    = Prob   .eq. FP

      Elastc    = lEmode .eq. 2
      call iload ( nb, 0, hEstat, 1 )

      !-----------------------------------------------------------------
      ! Find a feasible point.
      ! If the constraints are infeasible, minimize the sum of the
      ! elastic variables, subject to keeping the non-elastic variables
      ! feasible.  Elastic variables can move outside their bounds.
      !-----------------------------------------------------------------
      if (Elastc) then          ! make the nonelastics feasible
         lEmod0 = lEmode
         lvlIn0 = lvlInf
      else                      ! make everything feas or min sum
         lEmod0 = 1
         lvlIn0 = 2
      end if

      gotR    = .false.
      zngQP   = 0               ! No objective in phase 1
      znnH    = 0               ! No Hessian either

      iw(MnrHdP) = 1            ! Refresh print   heading.
      iw(MnrHdS) = 1            ! Refresh summary heading

      if (gotFac .eq. 0) then
         LUreq  = 1
      else
         LUreq  = 0
      end if
      typeLU  = BS
      LUok    = .true.
      done    = .false.

*     ==================================================================
*+    while (.not. done  .and.  LUok) do
  500 if    (.not. done  .and.  LUok) then

         needLU = LUreq .gt. 0

         call s2Bfac
     &      ( iExit, typeLU, needLU, newLU, newB,
     &        iObj, itn, PrtLvl, LUreq,
     &        m, mBS, n, nb, nnH, nS, nSwap,
     &        ne, nlocA, locA, indA, Acol,
     &        kBS, hs, bl, bu, blBS, buBS,
     &        nrhs0, nrhs, rhs, x, xBS,
     &        iy, iy1, y, y1, iw, leniw, rw, lenrw )
         if (iExit .ne. 0) go to 900

         needx  = needLU

         call s5QP
     &      ( inform, FP, ProbTag, Elastc, subopt,
     &        Hprod, Hprod1, QPlog, gotR, needLU, typeLU, needx,
     &        lenR, m, maxS, mBS, n, nb, nDegen, Hcalls,
     &        ngQP0, zngQP, ngObj0, ngObj, nnH0, znnH, nS,
     &        itQP, itQPmax, itQPtargt, itn, lEmod0, lvlIn0, PrtLvl,
     &        minimz, iObj, sclObj, ObjAdd, ObjFP, targtH, targtZ,
     &        tolFP, tolQP, tolx, nInf, sInf, wtInf, piNorm, rgNorm,
     &        ne, nlocA, locA, indA, Acol,
     &        hEtype, hEstat, hfeas, hs, kBS,
     &        Ascale, bl, bu, blBS, buBS,
     &        gBS, gObj, gQP, Hdx, pBS, pi, R, rc, rg,
     &        nrhs0, nrhs, rhs, lenx0, nx0, x0,
     &        x, xBS, xFreez,
     &        iy, iy1, y, y1, y2,
     &        cu, lencu, iu, leniu, ru, lenru,
     &        cw, lencw, iw, leniw, rw, lenrw )

         ! Check for trouble in s5QP.  Here are the possibilities:

         badLU   = inform .gt.  0 ! Fatal LU error
         found   = inform .eq.  0 ! QP solution found
         infeas  = inform .eq. -1 ! QP is infeasible
         unbnded = inform .eq. -2 ! QP is unbounded
         itnBig  = inform .eq. -3 ! Too many iterations

         if (badLU) then
            iExit = inform      ! Fatal LU error
            go to 900
         end if

         done = found .or. infeas .or. itnBig

         if (.not. done) then
            !===========================================================
            ! Trouble.
            ! Looks like the phase 1 problem is unbounded, which can
            ! occur only if a bad basis has given a large direction.
            !===========================================================
            ! Assume infeasible and factor with tighter tols.

            infeas = .true.

            write(str, 1000) itn
            call snPRNT( 3, str, iw, leniw )
            call s2tryLU
     &         ( itn, 22, nS, LUreq, LUok, typeLU,
     &           iw, leniw, rw, lenrw )
         end if
         go to 500
      end if
*+    end while
*     ==================================================================

      if (.not. found                        ) go to 800 ! Itns or infs
      if (nInf   .gt. 0  .and.  lEmode .ne. 1) go to 800 ! Infeas

      if (getOpt) then
!        call dcopy ( nb, x, 1, xFreez, 1 ) ! Save feasible nonelastics

!        if (PrtLvl .ge. 1  .and.  PrtLvl .lt. 10) then
         if (PrtLvl .ge. 10) then
            if ( Elastc ) then
               write(str, 1200) itn
               call snPRNT( 3, str, iw, leniw )
            else
               write(str, 1300) itn
               call snPRNT( 3, str, iw, leniw )
            end if
         end if
      end if

      !=================================================================
      ! Get optimal.
      !=================================================================
      if (Elastc) then
         ! Relax
      else if (lEmode .eq. 1) then ! Elastc mode never needed
         lvlInf = 1
      end if

      iw(MnrHdP) = 1               ! Refresh print heading.
      iw(MnrHdS) = 1

      indefQP = .false.

      LUreq   = 0
      typeLU  = BT
      LUok    = .true.
      done    = FPonly

*     ==================================================================
*+    while (.not. done  .and.  LUok) do
  600 if    (.not. done  .and.  LUok) then
         !--------------------------------------------------------------
         ! Compute and factorize the initial Z'HZ.
         ! The basis is refactorized if necessary.
         !--------------------------------------------------------------
         if (.not. gotR  .or.  LUreq .gt. 0) then
            call s5getR
     &         ( iExit, Hprod, Hprod1,
     &           Hcalls, gotR, typeLU, LUreq,
     &           itn, lenR, m, mBS, n, nb,
     &           nnH, nS, PrtLvl, minimz, iObj, targtH, targtZ,
     &           ne, nlocA, locA, indA, Acol,
     &           hs, kBS, bl, bu, blBS, buBS, R,
     &           nrhs0, nrhs, rhs,
     &           x, xBS, iy, iy1, y, y1, y2,
     &           cu, lencu, iu, leniu, ru, lenru,
     &           cw, lencw, iw, leniw, rw, lenrw )
            if (iExit .ne. 0) go to 900
         end if

         !--------------------------------------------------------------
         ! Solve the QP.
         !--------------------------------------------------------------
         LUreq  = 0

         call s5QP
     &      ( inform, Prob, ProbTag, Elastc, subopt,
     &        Hprod, Hprod1, QPlog, gotR, needLU, typeLU, needx,
     &        lenR, m, maxS, mBS, n, nb, nDegen, Hcalls,
     &        ngQP0, ngQP, ngObj0, ngObj, nnH0, nnH, nS,
     &        itQP, itQPmax, itQPtargt, itn, lEmode, lvlInf, PrtLvl,
     &        minimz, iObj, sclObj, ObjAdd, ObjQP, targtH, targtZ,
     &        tolFP, tolQP, tolx, nInf, sInf, wtInf, piNorm, rgNorm,
     &        ne, nlocA, locA, indA, Acol,
     &        hEtype, hEstat, hfeas, hs, kBS,
     &        Ascale, bl, bu, blBS, buBS,
     &        gBS, gObj, gQP, Hdx, pBS, pi, R, rc, rg,
     &        nrhs0, nrhs, rhs, lenx0, nx0, x0,
     &        x, xBS, xFreez,
     &        iy, iy1, y, y1, y2,
     &        cu, lencu, iu, leniu, ru, lenru,
     &        cw, lencw, iw, leniw, rw, lenrw )

         ! Check for trouble in s5QP.  Here are the possibilities:

         badLU   = inform .gt.  0 ! Fatal LU error
         found   = inform .eq.  0 ! QP solution found
         infeas  = inform .eq. -1 ! QP is infeasible
         unbnded = inform .eq. -2 ! QP is unbounded
         itnBig  = inform .eq. -3 ! Too many iterations
         weakMin = inform .eq. -4 ! Weak QP minimizer
         nSBig   = inform .eq. -5 ! Too many superbasics
         indefH  = inform .eq. -6 ! QP Hessian not positive semidefinite
         badZg   = inform .eq. -7 ! Z'g could not be made small enough
         badZ    = inform .eq. -8 ! Ill-conditioned Z

         if (badLU) then
            iExit = inform
            go to 900
         end if

         indefQP = indefH .and. indefQP

         done    = found  .or. unbnded .or. infeas .or. nSbig
     &                    .or. indefQP .or. itnBig .or. weakMin

         if (done) then
            !-----------------------------------------------------------
            ! Relax, we are finished.
            !-----------------------------------------------------------
         else
            !-----------------------------------------------------------
            ! Numerical trouble in s5QP.
            !-----------------------------------------------------------
            gotR  = .false.

            if (infeas) then
               !--------------------------------------------------------
               ! The nonelastics are infeasible. This should not happen.
               ! Phase 1 has already found a feasible point for the
               ! nonelastics, so the basis must be ill-conditioned.
               ! Refactorize with tighter tols and restart at the known
               ! feasible point.  Reduce the feasibility tol to try and
               ! prevent repeats.
               !-------------------------------------------------------
               write(str, 1100) itn
               call snPRNT( 3, str, iw, leniw )
               call s2tryLU
     &            ( itn, 22, nS, LUreq, LUok, typeLU,
     &              iw, leniw, rw, lenrw )

               elastc = .false.
!              call dcopy  ( nb, xFreez, 1, x, 1 )

            else if (indefH) then
               !--------------------------------------------------------
               ! Indefinite Z'HZ.  Could be an ill-conditioned Z'HZ.
               !--------------------------------------------------------
               write(str, 1600) itn
               call snPRNT( 3, str, iw, leniw )

               indefQP = indefH
               targtH  = one/(eps*eps)
               call s2tryLU
     &            ( itn, 26, nS, LUreq, LUok, typeLU,
     &              iw, leniw, rw, lenrw )

            else if (badZg) then
               !--------------------------------------------------------
               ! Large Z'g.        Probably an ill-conditioned Z.
               !--------------------------------------------------------
               write(str, 1700) itn
               call snPRNT( 23, str, iw, leniw )

               call s2tryLU
     &            ( itn, 21, nS, LUreq, LUok, typeLU,
     &              iw, leniw, rw, lenrw )

            else if (badZ) then
               !--------------------------------------------------------
               ! condZ > targtZ  while forming Z'HZ for a freq. check.
               ! Refactorize B, possibly with a reduced factor tol. If
               ! the factor tol is already tight, accept Z, however bad.
               !--------------------------------------------------------
               write(str, 1800) itn
               call snPRNT( 3, str, iw, leniw )

               call s2tryLU
     &            ( itn, 25, nS, LUreq, LUok, typeLU,
     &              iw, leniw, rw, lenrw )
               if (.not. LUok) then
                  targtZ = plInfy
                  LUok   = .true.
               end if
            end if
         end if

         go to 600
      end if
*+    end while
*     ==================================================================

  800 if (found) then
         if (nInf .eq. 0) then
            if (FPonly) then
               iExit =  2       ! feasible
            else
               iExit =  1       ! optimal
            end if
         else
            iExit = 14          ! optimal in elastic mode
         end if
      else if (infeas ) then
            iExit = 11          ! infeasible nonelastics
      else if (unbnded) then
            iExit = 21          ! unbounded
      else if (itnBig ) then
            iExit = 31          ! too many iterations
      else if (weakMin) then
            iExit =  4          ! weak minimizer
      else if (nSBig  ) then
            iExit = 33          ! too many superbasics
      else if (indefQP) then
            iExit = 53          ! Hessian indefinite
      else
            iExit = 41          ! Current point cannot be improved
      end if

  900 return

 1000 format(' Itn', i7, ': Infeasible nonelastics in feasibility',
     &                   ' phase')
 1100 format(' Itn', i7, ': Infeasible nonelastics in QP optimality',
     &                   ' phase')
 1200 format(' Itn', i7, ': Feasible non-elastics')
 1300 format(' Itn', i7, ': Feasible constraints')

 1600 format(' Itn', i7, ': Indefinite reduced QP Hessian --',
     &                   ' refactor the basis')
 1700 format(' Itn', i7, ': Large reduced gradient')
 1800 format(' Itn', i7, ': Ill-conditioned null space')

      end ! subroutine s5sQP

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s5sQN
     &   ( iExit, Hprod, Hprod1, QPlog, gotR, Prob,
     &     lenR, m, maxS, mBS, n, nb, nDegen, Hcalls,
     &     ngQP0, ngQP, ngObj0, ngObj, nnH0, nnH, nS,
     &     itQP, itQPmax, itn,
     &     minimz, iObj, sclObj, ObjAdd, ObjQP,
     &     tolFP, tolQP, tolx, nInf, sInf, wtInf, piNorm,
     &     ne, nlocA, locA, indA, Acol,
     &     hEtype, hEstat, hfeas, hs, kBS,
     &     Ascale, bl, bu, blBS, buBS,
     &     gBS, gObj, gQP, Hdx, pBS, pi, R, rc, rg, rg2,
     &     nrhs0, nrhs, rhs, lenx0, nx0, x0,
     &     x, xBS, xFreez,
     &     iy, iy1, y, y1, y2,
     &     cu, lencu, iu, leniu, ru, lenru,
     &     cw, lencw, iw, leniw, rw, lenrw )

      implicit
     &     none
      external
     &     Hprod, Hprod1, QPlog
      logical
     &     gotR
      integer
     &     Hcalls, iExit, iObj, itn, itQP, itQPmax, lenR, lenx0, lencu,
     &     leniu, lenru, lencw, leniw, lenrw, m, maxS, mBS, minimz,
     &     n, nb, ne, nInf, ngQP0, ngQP, nnH0, nnH, nS, nDegen,
     &     ngObj0, ngObj, nlocA, nrhs0, nrhs, nx0, Prob,
     &     hfeas(mBS), kBS(mBS), hEtype(nb), hEstat(nb), hs(nb),
     &     locA(nlocA), indA(ne), iy(nb), iy1(nb), iu(leniu), iw(leniw)
      double precision
     &     sclObj, ObjAdd, ObjQP, sInf, wtInf, piNorm, tolFP, tolQP,
     &     tolx, Acol(ne), Ascale(nb), bl(nb), bu(nb), blBS(mBS),
     &     buBS(mBS), gBS(mBS), gObj(ngObj0), gQP(ngQP0), Hdx(nnH0),
     &     pBS(nb), pi(m), R(lenR), rc(nb), rg(maxS), rg2(maxS),
     &     rhs(nrhs0), x(nb), x0(lenx0), xBS(mBS), xFreez(nb),
     &     y(nb), y1(nb), y2(nb), ru(lenru), rw(lenrw)
      character
     &     cu(lencu)*8, cw(lencw)*8

      !=================================================================
      ! s5sQN   solves the current problem.  An initial basis is assumed
      ! to be specified by nS, hs, x and the superbasic parts of kBS.
      ! In particular, there must be nS values hs(j) = 2, and the
      ! corresponding j's must be listed in kBS(m+1) thru kBS(m+nS).
      ! The ordering in kBS matches the reduced Hessian R (if any).
      !
      ! 01 Nov 2003: First version of s5sQN based on s5sQP.
      ! 08 Mar 2004: gotFac, gotHes implemented.
      ! 16 May 2006: Explicit target itQP added.
      ! 28 Apr 2011: Some cosmetic restructuring.
      !=================================================================
      character
     &     ProbTag*20, str*120
      logical
     &     badLU, badZg, badZ, done, Elast0, Elastc, found, FPonly,
     &     getOpt, indefH, infeas, itnBig, LUok, needLU, needx, newB,
     &     newLU, nSBig, unbnded, weakMin, ZitnBig
      integer
     &     gotFac, inform, itQPtargt, lEmod0, lEmode, LUreq,
     &     lvlIn0, lvlInf, MnrHdP, MnrHdS, MnrPrt, nSwap,
     &     PrtLvl, subopt, typeLU, zngQP
      double precision
     &     condHz, flmax, ObjFP, plInfy, targtZ, Zcndbd
      !-----------------------------------------------------------------
      integer            BS        , BT
      parameter         (BS     = 2, BT     = 3)
      integer            FP,         LP,         QP
      parameter         (FP     = 0, LP     = 1, QP = 2)
      integer            No
      parameter         (No     =-1)

      parameter         (MnrHdP = 223) ! >0 => Mnr heading for iPrint
      parameter         (MnrHdS = 224) ! >0 => Mnr heading for iSumm

      double precision   zero
      parameter         (zero   = 0.0d+0)
      !-----------------------------------------------------------------
      lEmode    = iw( 56) ! >0    => use elastic mode
      lvlInf    = iw( 73) ! Elastic option
      MnrPrt    = iw( 93) ! Minor print level
      gotFac    = iw(230) ! >0 => Save the LU factors

      flmax     = rw(  8) ! est. of the largest pos. real
      Zcndbd    = rw( 86) ! bound on the condition of Z

      PrtLvl    = MnrPrt
      subopt    = No

      plInfy    = flmax
      ObjQP     = zero
      targtZ    = Zcndbd
      condHz    = zero
      itQPtargt = itQPmax

      ProbTag   = 'linear constraints'

      getOpt    = Prob   .eq. LP   .or.  Prob .eq. QP
      FPonly    = Prob   .eq. FP

      Elastc    = lEmode .eq. 2
      call iload ( nb, 0, hEstat, 1 )

      !-----------------------------------------------------------------
      ! Find a feasible point.
      ! If the constraints are infeasible, minimize the sum of the
      ! elastic variables, subject to keeping the non-elastic variables
      ! feasible.  Elastic variables can move outside their bounds.
      !-----------------------------------------------------------------
      lvlIn0 = 2                ! local value of lvlInf
      lEmod0 = 1                ! local value of lEmode
      Elast0 = .false.          ! local value of Elastc

      zngQP  = 0                ! No objective in phase 1
      needx  = .false.          ! needLU will ask for new factors.

      if (gotFac .eq. 0) then
         LUreq  = 1
      else
         LUreq  = 0
      end if
      typeLU = BS
      LUok   = .true.
      done   = .false.

*     ==================================================================
*+    while (.not. done  .and.  LUok) do
  500 if    (.not. done  .and.  LUok) then

         needLU = LUreq .gt. 0

         call s2Bfac
     &      ( iExit, typeLU, needLU, newLU, newB,
     &        iObj, itn, PrtLvl, LUreq,
     &        m, mBS, n, nb, nnH, nS, nSwap,
     &        ne, nlocA, locA, indA, Acol,
     &        kBS, hs, bl, bu, blBS, buBS,
     &        nrhs0, nrhs, rhs, x, xBS,
     &        iy, iy1, y, y1, iw, leniw, rw, lenrw )
         if (iExit .ne. 0) go to 900

         needx = needLU

         call s5QN
     &      ( inform, FP, ProbTag, Elast0, subopt,
     &        Hprod, Hprod1, QPlog, gotR, needLU, typeLU, needx,
     &        lenR, m, maxS, mBS, n, nb, nDegen, Hcalls,
     &        ngQP0, zngQP, ngObj0, ngObj, nnH0, nnH, nS,
     &        itQP, itQPmax, itQPtargt, itn, lEmod0, lvlIn0, PrtLvl,
     &        minimz, iObj, sclObj, ObjAdd, ObjFP, condHz,
     &        targtZ, tolFP, tolQP, tolx, nInf, sInf, wtInf, piNorm,
     &        ne , nlocA, locA, indA, Acol,
     &        hEtype, hEstat, hfeas, hs, kBS,
     &        Ascale, bl, bu, blBS, buBS,
     &        gBS, gObj, gQP, Hdx, pBS, pi, R, rc, rg, rg2,
     &        nrhs0, nrhs, rhs, lenx0, nx0, x0,
     &        x, xBS, xFreez,
     &        iy, iy1, y, y1, y2,
     &        cu, lencu, iu, leniu, ru, lenru,
     &        cw, lencw, iw, leniw, rw, lenrw )

         ! Check for trouble in s5QN.  Here are the possibilities:

         badLU   = inform .gt.  0 ! Fatal LU error
         found   = inform .eq.  0 ! QP solution found
         infeas  = inform .eq. -1 ! QP is infeasible
         unbnded = inform .eq. -2 ! QP is unbounded
         itnBig  = inform .eq. -3 ! Too many iterations

         if (badLU) then
            iExit = inform      ! Fatal LU error
            go to 900
         end if

         done = found .or. infeas .or. itnBig

         if (.not. done) then
            !===========================================================
            ! Trouble.
            ! The phase 1 was unbounded, which can only occur if a bad
            ! basis has given a large search direction.
            !===========================================================
            ! Assume infeasible constraints and factor with tighter tols.

            infeas = .true.

            write(str, 1000) itn
            call snPRNT( 3, str, iw, leniw )
            call s2tryLU
     &         ( itn, 22, nS, LUreq, LUok, typeLU,
     &           iw, leniw, rw, lenrw )
         end if
         go to 500
      end if
*+    end while
*     ==================================================================

      if (.not. found) go to 800 ! Itns or infeasible

      if (getOpt) then
!        call dcopy ( nb, x, 1, xFreez, 1 ) ! Save feasible nonelastics

!        if (PrtLvl .ge. 1  .and.  PrtLvl .lt. 10) then
         if (PrtLvl .ge. 10) then
            if ( Elastc ) then
               write(str, 1200) itn
               call snPRNT( 3, str, iw, leniw )
            else
               write(str, 1300) itn
               call snPRNT( 3, str, iw, leniw )
            end if
         end if
      end if

      !=================================================================
      ! Get optimal.
      !=================================================================
      LUreq  = 0
      typeLU = BT
      LUok   = .true.
      done   = FPonly

*     ==================================================================
*+    while (.not. done  .and.  LUok) do
  600 if    (.not. done  .and.  LUok) then
         !--------------------------------------------------------------
         ! Refactorize the basis if necessary.
         !--------------------------------------------------------------
         needLU = LUreq .gt. 0

         if ( needLU ) then
            call s2Bfac
     &         ( iExit, typeLU, needLU, newLU, newB,
     &           iObj, itn, PrtLvl, LUreq,
     &           m, mBS, n, nb, nnH, nS, nSwap,
     &           ne, nlocA, locA, indA, Acol,
     &           kBS, hs, bl, bu, blBS, buBS,
     &           nrhs0, nrhs, rhs, x, xBS,
     &           iy, iy1, y, y1, iw, leniw, rw, lenrw )
            if (iExit .ne. 0) goto 900
         end if

         ! How do we update R if the superbasics change?
         !--------------------------------------------------------------
         ! Solve the QP using a quasi-Newton method.
         !--------------------------------------------------------------
         iw(MnrHdP) = 1         ! Refresh print heading.
         iw(MnrHdS) = 1

         LUreq  = 0

         if ( Elastc ) then
            lvlInf = 1 ! W1 = 1.0, W2 = wtInf  Elastic Phase 2 composite obj
         end if

        call s5QN
     &      ( inform, Prob, ProbTag, Elastc, subopt,
     &        Hprod, Hprod1, QPlog, gotR, needLU, typeLU, needx,
     &        lenR, m, maxS, mBS, n, nb, nDegen, Hcalls,
     &        ngQP0, ngQP, ngObj0, ngObj, nnH0, nnH, nS,
     &        itQP, itQPmax, itQPtargt, itn, lEmode, lvlInf, PrtLvl,
     &        minimz, iObj, sclObj, ObjAdd, ObjQP, condHz,
     &        targtZ, tolFP, tolQP, tolx, nInf, sInf, wtInf, piNorm,
     &        ne , nlocA, locA, indA, Acol,
     &        hEtype, hEstat, hfeas, hs, kBS,
     &        Ascale, bl, bu, blBS, buBS,
     &        gBS, gObj, gQP, Hdx, pBS, pi, R, rc, rg, rg2,
     &        nrhs0, nrhs, rhs, lenx0, nx0, x0,
     &        x, xBS, xFreez,
     &        iy, iy1, y, y1, y2,
     &        cu, lencu, iu, leniu, ru, lenru,
     &        cw, lencw, iw, leniw, rw, lenrw )

         ! Check for trouble in s5QN.  Here are the possibilities:

         badLU   = inform .gt.  0 ! Fatal LU error
         found   = inform .eq.  0 ! QP solution found
         infeas  = inform .eq. -1 ! QP is infeasible
         unbnded = inform .eq. -2 ! QP is unbounded
         itnBig  = inform .eq. -3 ! Too many iterations
         weakMin = inform .eq. -4 ! Void, should not happen
         nSBig   = inform .eq. -5 ! Too many superbasics
         indefH  = inform .eq. -6 ! Void, should not happen
         badZg   = inform .eq. -7 ! Void, should not happen
         badZ    = inform .eq. -8 ! Ill-conditioned Z
         ZitnBig = inform .eq. -9 ! Too many subspace iterations

         if (badLU) then
            iExit = inform
            go to 900
         end if

         done    = found  .or. unbnded .or. infeas .or. nSbig .or.
     &             itnBig .or. weakMin .or. indefH .or. badZg .or.
     &             ZitnBig

         if (done) then
            !-----------------------------------------------------------
            ! Relax, we are finished.
            !-----------------------------------------------------------
         else
            !-----------------------------------------------------------
            ! Numerical trouble in s5QN.
            !-----------------------------------------------------------
            gotR  = .false.

            if (infeas) then
               !--------------------------------------------------------
               ! The nonelastics are infeasible. This should not happen.
               ! Phase 1 has already found a feasible point for the
               ! nonelastics, so the basis must be ill-conditioned.
               ! Refactorize with tighter tols and restart at the known
               ! feasible point.  Reduce the feasibility tol to try and
               ! prevent repeats.
               !--------------------------------------------------------
               write(str, 1100) itn
               call snPRNT( 3, str, iw, leniw )
               call s2tryLU
     &            ( itn, 22, nS, LUreq, LUok, typeLU,
     &              iw, leniw, rw, lenrw )

               elastc = .false.
!              call dcopy  ( nb, xFreez, 1, x, 1 )

            else if (badZ) then
               !--------------------------------------------------------
               ! Numerical error.  Probably a bad basis.
               ! Refactorize B, possibly with a reduced factor tol. If
               ! the factor tol is already tight, accept Z, however bad.
               !--------------------------------------------------------
               write(str, 1800) itn
               call snPRNT( 3, str, iw, leniw )
               call s2tryLU
     &            ( itn, 25, nS, LUreq, LUok, typeLU,
     &              iw, leniw, rw, lenrw )
               if (.not. LUok) then
                  targtZ = plInfy
                  LUok   = .true.
               end if
            end if
         end if

         go to 600
      end if
*+    end while
*     ==================================================================

  800 if (found) then
         if (nInf .eq. 0) then
            iExit =  1          ! optimal
         else
            iExit = 14          ! optimal
         end if
      else if (infeas ) then
            iExit = 11          ! infeasible nonelastics
      else if (unbnded) then
            iExit = 21          ! unbounded
      else if (itnBig ) then
            iExit = 31          ! too many iterations
      else if (nSBig) then
            iExit = 33          ! too many superbasics
      else
            iExit = 41          ! Current point cannot be improved
      end if

  900 return

 1000 format(' Itn', i7, ': Infeasible nonelastics in feasibility',
     &                   ' phase')
 1100 format(' Itn', i7, ': Infeasible nonelastics in QN optimality',
     &                   ' phase')
 1200 format(' Itn', i7, ': Feasible non-elastics')
 1300 format(' Itn', i7, ': Feasible constraints')

 1600 format(' Itn', i7, ': Indefinite reduced Hessian')
 1700 format(' Itn', i7, ': Large reduced gradient')
 1800 format(' Itn', i7, ': Ill-conditioned null space')

      end ! subroutine s5sQN

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s5Stat
     &   ( Status, iw, leniw)

      implicit
     &     none

      integer
     &     Status, leniw, iw(leniw)

      !=================================================================
      ! s5Stat fetches the call-status for the sqOpt user-defined
      ! matrix-vector product.
      !
      ! 16 Jun 2008: First version of s5Stat.
      !=================================================================
      character
     &     str*80
      integer
     &     qpStat
      !-----------------------------------------------------------------
      parameter         (qpStat = 235) ! QP user-routine call-status
      !-----------------------------------------------------------------

      if (     iw(qpStat) .eq. 0) then
         ! Standard call

         Status     =  0
      else if (iw(qpStat) .lt. 0) then
         ! First call

         Status     =  1
         iw(qpStat) =  0
      else if (iw(qpStat) .ge. 2) then
         ! Last orders please

         Status     = iw(qpStat)
         iw(qpStat) = -1
      else
         Status     = iw(qpStat)
         write(str, 9999) Status
         call snPRNT( 3, str, iw, leniw )
      end if

      return

 9999 format(' XXX  user-function call-status not recognized.',
     &       ' Requested status =', i6)

      end ! subroutine s5Stat
