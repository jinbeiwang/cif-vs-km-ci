/*=======================================================================================
  CIF PERCENTILE 95% CI BY INVERTING THE CIF CONFIDENCE LIMITS -- WORKED EXAMPLE
  =======================================================================================
  PURPOSE
    A self-contained, hand-checkable example that starts from an ALREADY-COMPUTED CIF
    dataset (a list of time / CIF / CIF_LCL / CIF_UCL rows) and produces the percentile
    point estimate and its 95% CI by inverting the confidence limits. Every step is
    printed so you can trace the arithmetic by hand.

    This is the same mechanism as SAS's own KM quartiles (Brookmeyer-Crowley inversion),
    applied to a NON-DECREASING step function. Because F is increasing while S is
    decreasing, the two confidence limits trade roles compared to the KM branch:
        lower bound q_L = first t where the UPPER limit reaches p
        upper bound q_U = first t where the LOWER limit reaches p

  INPUT
    The CIF dataset below (work.cif_demo) is a small, FIXED table. In production this
    is exactly what OUTCIF= returns from PROC LIFETEST (columns: time, CIF, CIF_LCL,
    CIF_UCL). The numbers are chosen so the crossing is easy to see by eye.

  OUTPUT
    Table A  the input CIF step function, one row per event time
    Table B  the run of F, F_lcl, F_ucl with cumulative maxima (plateaus) shown
    Table C  the final percentile table (p=0.25/0.50/0.75): q_hat, q_lcl, q_ucl,
             q_prev, F_at_q, and the two self checks chk_order / chk_hit

  Version: SAS 9.4M7 / SAS/STAT 15.2. ASCII English only (server-safe).
=======================================================================================*/

options nodate nonumber linesize=160 formdlim='-';

/*==== 1. INPUT: a fixed CIF step-function dataset ===================================*/
/*  Columns required by the inversion macro (names matter):
        pop      group id (here a single group "Demo")
        time     event time at which the step jumps (right-continuous convention)
        F        CIF estimate                F(t)
        F_lcl    lower 95% pointwise limit   CIF_LCL(t)
        F_ucl    upper 95% pointwise limit   CIF_UCL(t)
    The rows are sorted by pop time. There is NO t=0 row with a missing F -- that row is
    dropped upstream (OUTCIF= carries it; see cif_km_quantile_demo.sas).                 */

data cif_demo;
  input time F F_lcl F_ucl;
  length pop $8;
  pop = 'Demo';
  label time  = 'Event time (months)'
        F     = 'CIF(t)'
        F_lcl = 'CIF lower 95% limit'
        F_ucl = 'CIF upper 95% limit';
  datalines;
1.0  0.08  0.02  0.22
2.0  0.15  0.06  0.34
3.0  0.21  0.10  0.41
4.0  0.27  0.14  0.48
5.0  0.32  0.18  0.55
6.0  0.36  0.21  0.60
7.0  0.39  0.24  0.64
8.0  0.41  0.26  0.67
;
run;

title "A. Input CIF step function (work.cif_demo)";
proc print data=cif_demo noobs label; run;

/*==== 2. THE INVERSION ==============================================================*/
/*  Parameters: the percentiles to invert and the confidence level (alpha is informational
    here because the limits are already computed; it documents the 95% setting).        */
%let pslist = 0.25,0.50,0.75;
%let nq     = 3;

/*  INVERSION DEFINITION (for a NON-DECREASING F):
        q_hat = inf{ t : F(t)     >= p }     first time where the CIF reaches p
        q_lcl = inf{ t : F_UCL(t) >= p }     first time where the UPPER limit reaches p
        q_ucl = inf{ t : F_LCL(t) >= p }     first time where the LOWER limit reaches p

    Because F, F_lcl and F_ucl are non-decreasing step functions, each expression is
    simply "scan the rows top-down and take the FIRST time that satisfies the condition".
    No root finding, no interpolation. The scan is done with retained temporary arrays,
    reset per group (BY pop).                                                           */

data cif_inverted;
  set cif_demo;
  by pop;

  array pctl[&nq.] _temporary_ (&pslist.);
  array tq[&nq.]   _temporary_;       /* q_hat : first t with F     >= p */
  array tl[&nq.]   _temporary_;       /* q_lcl : first t with F_ucl >= p */
  array tu[&nq.]   _temporary_;       /* q_ucl : first t with F_lcl >= p */
  array tp[&nq.]   _temporary_;       /* time just before the crossing    */
  array lv[&nq.]   _temporary_;       /* F(q_hat), the level at crossing  */

  retain plateau plateau_lcl plateau_ucl maxtime;

  if first.pop then do;
    do i = 1 to &nq.;
      tq[i] = .; tl[i] = .; tu[i] = .; tp[i] = .; lv[i] = .;
    end;
    plateau = .; plateau_lcl = .; plateau_ucl = .; maxtime = .;
  end;

  maxtime     = max(maxtime,     time);
  plateau     = max(plateau,     F);
  plateau_lcl = max(plateau_lcl, F_lcl);
  plateau_ucl = max(plateau_ucl, F_ucl);

  do i = 1 to &nq.;
    if tq[i] = . then do;
      if F >= pctl[i] then do; tq[i] = time; lv[i] = F; end;
      else tp[i] = time;               /* still below p -> remember this time */
    end;
    if tl[i] = . and F_ucl >= pctl[i] then tl[i] = time;
    if tu[i] = . and F_lcl >= pctl[i] then tu[i] = time;
  end;

  if last.pop then do i = 1 to &nq.;
    percentile = pctl[i];
    q_hat      = tq[i];
    q_lcl      = tl[i];
    q_ucl      = tu[i];
    q_prev     = tp[i];
    F_at_q     = lv[i];

    /* Self check 1: when all three exist, q_lcl <= q_hat <= q_ucl must hold.  */
    chk_order = .;
    if n(q_hat, q_lcl, q_ucl) = 3 then chk_order = (q_lcl <= q_hat <= q_ucl);

    /* Self check 2: the level at the crossing is >= p (right-continuity).     */
    chk_hit = .;
    if F_at_q ne . then chk_hit = (F_at_q >= percentile);

    output;
  end;

  keep pop percentile q_hat q_lcl q_ucl q_prev F_at_q
       plateau plateau_lcl plateau_ucl maxtime chk_order chk_hit;
  format percentile 4.2 F_at_q 5.3;
run;

/*==== 3. RESULTS ====================================================================*/

title "B. Cumulative maxima (plateaus) used to flag non-estimable cases";
proc sql;
  create table b_show as
  select time, F, F_lcl, F_ucl,
         max(F)     over (order by time) as cumF,
         max(F_lcl) over (order by time) as cumL,
         max(F_ucl) over (order by time) as cumU
  from cif_demo;
quit;
proc print data=b_show noobs; run;

title "C. Percentile inversion result (months)";
proc print data=cif_inverted noobs label;
  var pop percentile q_hat q_lcl q_ucl q_prev F_at_q
      plateau plateau_lcl plateau_ucl maxtime chk_order chk_hit;
run;

/*  Non-estimable handling: if q_hat is missing, the CIF plateau never reached p (report
    NE). If q_hat exists but q_lcl or q_ucl is missing, that bound is non-estimable: the
    corresponding confidence limit plateau never reached p (report NE for that bound and
    footnote the plateau value and max follow-up). The columns plateau / plateau_lcl /
    plateau_ucl / maxtime are carried for exactly that documentation.

    EXPECTED RESULT for this dataset (months), all hand-checkable:

      p=0.25:  q_hat = 4  (F(4)=0.27 is the first F >= 0.25)
               q_lcl = 2  (F_ucl(2)=0.34 is the first F_ucl >= 0.25)
               q_ucl = 8  (F_lcl(8)=0.26 is the first F_lcl >= 0.25)
               => 2 <= 4 <= 8, chk_order = 1; F(4)=0.27 >= 0.25, chk_hit = 1

      p=0.50:  q_hat = .  (F plateau = 0.41 < 0.50  -> NE)
               q_lcl = 5  (F_ucl(5)=0.55 is the first F_ucl >= 0.50)
               q_ucl = .  (F_lcl plateau = 0.26 < 0.50 -> NE)

      p=0.75:  all three NE (all three plateaus < 0.75)

    NOTE the two self checks: chk_order = 1 whenever all three bounds exist (they must
    satisfy q_lcl <= q_hat <= q_ucl), and chk_hit = 1 whenever q_hat exists (the CIF at
    the crossing must be >= p). Both must be 1 in every row where they are computed.     */

title;
/*=======================================================================================
  HOW THIS MAPS TO THE PRODUCTION PROGRAM
    - cif_demo above <-> the OUTCIF= dataset of PROC LIFETEST (after normalising column
      names and dropping the t=0 row with a missing CIF). See cif_km_quantile_demo.sas
      macro mkFailFromCif for that normalisation.
    - The data step "cif_inverted" <-> the invCurve macro in cif_km_quantile_demo.sas
      (generalised there to arbitrary pop groups and any list of percentiles).
=======================================================================================*/
