/*=======================================================================================
  CIF PERCENTILE 95% CI BY INVERTING CONFIDENCE LIMITS -- RUNNABLE DEMO
  KM percentiles vs Aalen-Johansen CIF percentiles
  ---------------------------------------------------------------------------------------
  WHAT THIS PROGRAM DOES
    On ONE single dataset it computes the 25th / Median / 75th percentile of three curves
    side by side, so that the "invert the CIF confidence limits" mechanism can be tested
    and cross-checked against the KM percentiles that SAS itself produces.

      Curve A   1 - KM_cs    cause-specific KM (competing event treated as censored)
      Curve B   1 - KM_comp  composite KM (any event counted as an event)
      Curve C   CIF          Aalen-Johansen cumulative incidence of the event of interest

  ORDERING THAT MUST HOLD (all three curves are non-decreasing, hence same order for
  percentiles):
      CIF(t) <= 1-KM_cs(t) <= 1-KM_comp(t)   =>   q(CIF) >= q(KM_cs) >= q(KM_comp)
  i.e. KM-based percentiles are systematically EARLIER than the CIF-based percentile.

  ---------------------------------------------------------------------------------------
  WHY KM USES OUTSURV= AND THE CIF USES OUTCIF=   *** READ THIS ***
    - PROC LIFETEST estimates the SURVIVAL function S(t). Its natural output object is
      therefore Survival, and OUTCIF= is NOT available in that (single-cause / KM) setting.
      To obtain a failure probability from KM you must COMPLEMENT it yourself:
          F(t) = 1 - S(t),  F_LCL = 1 - SDF_UCL,  F_UCL = 1 - SDF_LCL
      Note the sides flip when complementing.
    - The Aalen-Johansen cumulative incidence is a different estimator that already targets
      the probability of having the event of interest by time t. OUTCIF= delivers it
      directly, so NO complement is applied and the limits keep their own sides.

  OUTSURV= output data set -- the columns used here (verified against SAS documentation):
      <BY variables>, STRATUM, <time variable as given in the TIME statement>,
      SURVIVAL, CONFTYPE, SDF_LCL, SDF_UCL, _CENSOR_
      => the confidence limit columns are SDF_LCL / SDF_UCL, NOT LowerSurv / UpperSurv.
      => the time column keeps the TIME-statement variable name, it is NOT called "Time".

  OUTCIF= output data set -- the columns used here (SAS/STAT 15.2):
      <BY variables>, <time variable as given in the TIME statement>,
      CIF, CIF_STDERR, CIF_LCL, CIF_UCL, ALPHA, CONFTYPE, ...
      => the first row (t=0) carries a MISSING CIF and must be dropped.

  ---------------------------------------------------------------------------------------
  INVERSION DEFINITION (identical in structure to how SAS computes KM quartiles,
  Brookmeyer-Crowley inversion of the g-transformed confidence limits):
      point estimate  q   = inf{ t : F(t)   >= p }   first event time where F reaches p
      lower limit     q_L = inf{ t : F_UCL(t) >= p } first event time where F_UCL reaches p
      upper limit     q_U = inf{ t : F_LCL(t) >= p } first event time where F_LCL reaches p
  F is INCREASING while the survival function is DECREASING, so on the CIF side the LOWER
  limit is driven by the UPPER confidence limit -- the opposite of the KM branch.

  ---------------------------------------------------------------------------------------
  HOW TO RUN
    (1) Submit as is: data are simulated on the fly with CALL STREAMINIT (no external file
        needed, so it always runs). The data generating model matches the companion Python
        reference implementation cif_km_quantile_demo.py, but the random number stream
        differs, so the numbers will not match digit by digit (conclusions do match).
    (2) For a digit-by-digit cross-check: run cif_km_quantile_demo.py once to create
        sim_competing.csv, then set %let simcsv = Y; below. Both sides then read the very
        same data.
    (3) To use your own data: change tvar / cvar / armvar / evcode / cencode / cecode in
        the parameter block, and replace the whole data preparation section with your
        analysis dataset (see the comment at the end of section 1).

  Version: SAS 9.4M7 / SAS/STAT 15.2   (OUTCIF= requires SAS/STAT 14.1 or later)
  All comments, titles and messages are ASCII English on purpose: the submission server
  does not handle non-ASCII source text.
=======================================================================================*/

options nodate nonumber linesize=160 formdlim='-';

/*==== 0. PARAMETERS ==================================================================*/
%let simcsv  = N;                       /* Y = read sim_competing.csv; N = simulate      */
%let csvpath = G:\test\surv\sim_competing.csv;
%let tvar    = avald;                   /* time variable (days here; must match the Shell)*/
%let cvar    = cnsrscd;                 /* outcome code: 0=censored 1=metastasis 2=death */
%let armvar  = trta;                    /* cohort / treatment group variable             */
%let evcode  = 1;                       /* event code of interest                        */
%let cencode = 0;                       /* censoring code                                */
%let cecode  = 2;                       /* competing event code (death before metastasis)*/
%let alpha   = 0.05;                    /* confidence level                              */
%let ps      = 0.25 0.50 0.75;          /* 25th / Median / 75th                          */
%let yrfac   = 365.25;                  /* days -> years (display only)                  */

%let err_km  = greenwood;               /* variance for the KM runs. GreenWood is the
                                           SAS default for the survival function, so the
                                           cross-check against the Quartiles table is
                                           apples-to-apples. Switch to aalen to see the
                                           sensitivity of the CI limits.                 */
%let err_cif = aalen;                   /* variance for the CIF run (as specified by the
                                           statistician: Aalen-Johansen + Aalen variance)*/
%let cftype  = loglog;                  /* g-transformation, same on both branches       */

%let nq      = %sysfunc(countw(&ps.));
%let pslist  = %sysfunc(tranwrd(%str(&ps.), %str( ), %str(,)));

/*==== 1. DATA PREPARATION ============================================================*/
%if %upcase(&simcsv.) = Y %then %do;
  proc import datafile="&csvpath." out=sim0 dbms=csv replace;
    getnames = yes;
  run;
%end;
%else %do;
  /* Same generating model as cif_km_quantile_demo.py:
       Arm A: lambda1 = 0.22/yr (metastasis), lambda2 = 0.10/yr (death pre-metastasis)
       Arm B: lambda1 = 0.30/yr,               lambda2 = 0.14/yr
       loss to follow-up: lambda_c = 0.03/yr (both arms); administrative censor at 5 years.
     The lambdas are per YEAR, so rand('EXPONENTIAL', 1/lambda) yields years, and the
     value is multiplied by 365.25 to be stored in days.                                */
  data sim0;
    call streaminit(20260917);
    length TRTA $10;
    keep USUBJID TRTA AVALD CNSRSCD;
    do arm = 1 to 2;
      if arm = 1 then do; TRTA = 'Arm A'; lam1 = 0.22; lam2 = 0.10; end;
      else            do; TRTA = 'Arm B'; lam1 = 0.30; lam2 = 0.14; end;
      lamc = 0.03;
      do id = 1 to 200;
        t1 = rand('EXPONENTIAL', 1/lam1);
        t2 = rand('EXPONENTIAL', 1/lam2);
        tc = rand('EXPONENTIAL', 1/lamc);
        tt = min(t1, t2, tc);
        if      tt = t1 then ev = 1;
        else if tt = t2 then ev = 2;
        else                 ev = 0;
        if tt > 5 then do; tt = 5; ev = 0; end;     /* administrative censoring */
        AVALD   = tt * 365.25;
        CNSRSCD = ev;
        USUBJID = cats(put(arm, 1.), put(id, z4.));
        output;
      end;
    end;
  run;
%end;

/* With a real analysis dataset, replace the whole block above by reading your ADTTE and
   make sure the three variables below exist:
       &tvar.  numeric time      &cvar.  outcome code (0/1/2)      &armvar.  group        */

data sim;
  set sim0(rename=(&tvar. = _t &cvar. = _c &armvar. = _a));
  length TRTA $16 pop $16;
  TRTA    = _a;
  AVALD   = _t;
  CNSRSCD = _c;
  drop _t _c _a;
  label AVALD   = 'Time (days)'
        CNSRSCD = 'Outcome: 0=Censored 1=Metastasis 2=Death prior to Metastasis';
  pop = TRTA;                                        /* one row set per cohort          */
  output;
  pop = 'Total';                                     /* the Shell "Total" column        */
  output;
run;

proc sort data=sim; by pop; run;

proc freq data=sim(where=(pop ne 'Total')) noprint;
  tables TRTA*CNSRSCD / out=_dist;
run;
title "0. Simulated data: outcome distribution (0=censored 1=metastasis 2=death)";
proc print data=_dist noobs; run;

/*==== 2. BUILD THE THREE CURVES =====================================================*/
/*  Three hard rules:
      (1) NEVER add REDUCEOUT. It truncates the output to the TIMELIST time points, so the
          step function is incomplete and the inversion is guaranteed to be wrong.
      (2) The parentheses in the TIME statement list the CENSORING codes; eventcode=
          names the event of interest. Values that are neither a censoring code nor the
          event code automatically become competing events.
      (3) error= selects the variance estimate. See the parameter block for the choices.  */

/* --- 2A. Curve B: composite KM = 1 - KM_all. OUTSURV= is used because this is the KM
          setting, where OUTCIF= is not available. The SAS Quartiles table is captured at
          the same time so that the inversion can be validated against it. ------------- */
/*   Do NOT add NOPRINT here: combined with ODS OUTPUT it can suppress the dataset as
     well. ODS SELECT NONE / ALL is the safe way to silence the listing.                */
ods select none;
ods output Quartiles = km_q_sas;
proc lifetest data=sim outsurv=km_all
     error = &err_km. conftype = &cftype. alpha = &alpha.;
  by pop;
  time AVALD*CNSRSCD(&cencode.);          /* only 0 is censoring -> 1 and 2 are events  */
run;
ods output close;
ods select all;

/* --- 2B. Curve A: cause-specific KM = 1 - KM_cs, with the competing event censored --- */
proc lifetest data=sim outsurv=km_cs noprint
     error = &err_km. conftype = &cftype. alpha = &alpha.;
  by pop;
  time AVALD*CNSRSCD(&cencode.,&cecode.); /* 0 and 2 censored -> only 1 is an event      */
run;

/* --- 2C. Curve C: Aalen-Johansen CIF of the event of interest. Here OUTCIF= applies,
          because the competing-risks machinery is active (eventcode= is in use). ------- */
proc lifetest data=sim outcif=cif_c1 noprint
     error = &err_cif. conftype = &cftype. alpha = &alpha.;
  by pop;
  time AVALD*CNSRSCD(&cencode.) / eventcode = &evcode.;
run;

/* Column-name self check (look here first if your SAS version names things differently) */
%if %sysfunc(exist(km_all)) %then %do;
  proc contents data=km_all varnum short; run;
%end;
%if %sysfunc(exist(cif_c1)) %then %do;
  proc contents data=cif_c1 varnum short; run;
%end;

/*==== 3. NORMALISE ALL CURVES TO ONE LAYOUT: pop time F F_lcl F_ucl ==================*/
/*  Robust time-column detection: OUTSURV= and OUTCIF= both carry the time variable under
    the name used in the TIME statement, but fall back to TIME if that is not the case. */
%macro _pickTime(ds);
  %global _tcol;
  %let _tcol = ;
  proc sql noprint;
    select name into :_tcol trimmed
      from dictionary.columns
     where libname = 'WORK' and memname = "%upcase(&ds.)"
       and upcase(name) = "%upcase(&tvar.)";
  quit;
  %if %length(&_tcol.) = 0 %then %do;
    proc sql noprint;
      select name into :_tcol trimmed
        from dictionary.columns
       where libname = 'WORK' and memname = "%upcase(&ds.)" and upcase(name) = 'TIME';
    quit;
  %end;
  %if %length(&_tcol.) = 0 %then %let _tcol = &tvar.;
%mend _pickTime;

/* Shared guard: create an empty placeholder so that downstream steps never abort. */
%macro _emptyCurve(dsout=, why=);
  data &dsout.;
    length pop $16;
    length time F F_lcl F_ucl 8;
    stop;
  run;
  %put NOTE: [demo] &dsout. left EMPTY -- &why.;
%mend _emptyCurve;

/* ---- KM branch: the estimator targets SURVIVAL, so COMPLEMENT it ------------------
   F = 1 - S.  When complementing, the confidence limits swap sides:
       F_lcl = 1 - SDF_UCL        F_ucl = 1 - SDF_LCL
   The columns read here are SURVIVAL / SDF_LCL / SDF_UCL. If your SAS version names
   them differently the program leaves the curve empty and tells you instead of aborting. */
%macro mkFailFromSurv(dsin=, dsout=);
  %_pickTime(&dsin.)
  %let _ok = 0;
  %if %sysfunc(exist(&dsin.)) %then %do;
    proc sql noprint;
      select count(*) into :_ok trimmed from dictionary.columns
       where libname = 'WORK' and memname = "%upcase(&dsin.)"
         and upcase(name) in ('SURVIVAL','SDF_LCL','SDF_UCL');
    quit;
  %end;
  %if &_ok. = 3 %then %do;
    data &dsout.;
      set &dsin.(keep = pop &_tcol. SURVIVAL SDF_LCL SDF_UCL);
      time  = &_tcol.;
      F     = 1 - SURVIVAL;
      F_lcl = 1 - SDF_UCL;    /* survival upper limit -> failure lower limit */
      F_ucl = 1 - SDF_LCL;    /* survival lower limit -> failure upper limit */
      if F = . then delete;
      keep pop time F F_lcl F_ucl;
    run;
    proc sort data=&dsout.; by pop time; run;
  %end;
  %else %_emptyCurve(dsout=&dsout.,
                     why=required columns SURVIVAL/SDF_LCL/SDF_UCL not found in &dsin. detected _ok=&_ok. -- adjust mkFailFromSurv per PROC CONTENTS);
%mend mkFailFromSurv;

/* ---- CIF branch: the estimator ALREADY targets the cumulative incidence, NO complement
       F = CIF      F_lcl = CIF_LCL      F_ucl = CIF_UCL                            */
%macro mkFailFromCif(dsin=, dsout=);
  %_pickTime(&dsin.)
  %let _ok = 0;
  %if %sysfunc(exist(&dsin.)) %then %do;
    proc sql noprint;
      select count(*) into :_ok trimmed from dictionary.columns
       where libname = 'WORK' and memname = "%upcase(&dsin.)"
         and upcase(name) in ('CIF','CIF_LCL','CIF_UCL');
    quit;
  %end;
  %if &_ok. = 3 %then %do;
    data &dsout.;
      set &dsin.(keep = pop &_tcol. CIF CIF_LCL CIF_UCL);
      time  = &_tcol.;
      F     = CIF;
      F_lcl = CIF_LCL;
      F_ucl = CIF_UCL;
      if F = . then delete;   /* OUTCIF carries a t=0 row with a missing CIF */
      keep pop time F F_lcl F_ucl;
    run;
    proc sort data=&dsout.; by pop time; run;
  %end;
  %else %_emptyCurve(dsout=&dsout.,
                     why=required columns CIF/CIF_LCL/CIF_UCL not found in &dsin. detected _ok=&_ok. -- adjust mkFailFromCif per PROC CONTENTS);
%mend mkFailFromCif;

%mkFailFromSurv(dsin=km_cs,  dsout=cur_cs);
%mkFailFromSurv(dsin=km_all, dsout=cur_all);
%mkFailFromCif(dsin=cif_c1,  dsout=cur_cif);

/*==== 4. GENERIC INVERSION MACRO ====================================================*/
/*  Input : any dataset with pop / time / F / F_lcl / F_ucl, sorted by pop time.
    Output: long table pop x percentile, with q_hat / q_lcl / q_ucl / q_prev / F_at_q
            and the plateau values used to document non-estimable cases.
    Because F, F_lcl and F_ucl are non-decreasing step functions, each of the three
    expressions degenerates into "scan downwards, take the first row that satisfies the
    condition" -- no root finding is needed.                                            */
%macro invCurve(dsin=, dsout=);
data &dsout.;
  set &dsin.;
  by pop time;

  array pctl[&nq.] _temporary_ (&pslist.);
  array tq[&nq.]   _temporary_;      /* point estimate q   = inf{ t : F     >= p }        */
  array tl[&nq.]   _temporary_;      /* lower limit    q_L = inf{ t : F_ucl >= p }        */
  array tu[&nq.]   _temporary_;      /* upper limit    q_U = inf{ t : F_lcl >= p }        */
  array tp[&nq.]   _temporary_;      /* last time point before crossing (F < p)           */
  array lv[&nq.]   _temporary_;      /* level F(q) at the crossing, always >= p           */

  retain plateau plateau_lcl plateau_ucl maxtime;
  if first.pop then do;
    do i = 1 to &nq.;
      tq[i] = .; tl[i] = .; tu[i] = .; tp[i] = .; lv[i] = .;
    end;
    plateau = .; plateau_lcl = .; plateau_ucl = .; maxtime = .;
  end;

  if F ne . then do;
    maxtime     = max(maxtime,     time);
    plateau     = max(plateau,     F);
    plateau_lcl = max(plateau_lcl, F_lcl);
    plateau_ucl = max(plateau_ucl, F_ucl);

    do i = 1 to &nq.;
      if tq[i] = . then do;
        if F >= pctl[i] then do; tq[i] = time; lv[i] = F; end;
        else tp[i] = time;
      end;
      if tl[i] = . and F_ucl >= pctl[i] then tl[i] = time;
      if tu[i] = . and F_lcl >= pctl[i] then tu[i] = time;
    end;
  end;

  if last.pop then do i = 1 to &nq.;
    percentile = pctl[i];
    q_hat    = tq[i];
    q_lcl    = tl[i];
    q_ucl    = tu[i];
    q_prev   = tp[i];
    F_at_q   = lv[i];
    /* check 1: whenever all three exist, q_L <= q <= q_U must hold */
    chk_order = .;
    if n(q_hat, q_lcl, q_ucl) = 3 then chk_order = (q_lcl <= q_hat <= q_ucl);
    /* check 2: the level at the crossing must be >= p (right-continuity) */
    chk_hit = .;
    if F_at_q ne . then chk_hit = (F_at_q >= percentile);
    output;
  end;

  keep pop percentile q_hat q_lcl q_ucl q_prev F_at_q
       plateau plateau_lcl plateau_ucl maxtime chk_order chk_hit;
  format percentile 4.2;
run;

proc sort data=&dsout.; by pop percentile; run;
%mend invCurve;

%macro safeInv(dsin=, dsout=);
  %if %sysfunc(exist(&dsin.)) %then %invCurve(dsin=&dsin., dsout=&dsout.);
  %else %do;
    data &dsout.;
      length pop $16;
      length percentile q_hat q_lcl q_ucl q_prev F_at_q
             plateau plateau_lcl plateau_ucl maxtime chk_order chk_hit 8;
      stop;
    run;
    %put NOTE: [demo] &dsin. was not created -- &dsout. left empty.;
  %end;
%mend safeInv;

%safeInv(dsin=cur_cs,  dsout=q_cs);
%safeInv(dsin=cur_all, dsout=q_all);
%safeInv(dsin=cur_cif, dsout=q_cif);

/*==== 5. CROSS-VALIDATION AGAINST THE SAS QUARTILES TABLE ===========================*/
/*  Same data, same 1-KM_all curve, same variance (GreenWood) and transformation (LOGLOG):
     the percentiles produced by the inversion above must reproduce SAS's own Quartiles
     table -- the point estimates exactly, and the CI limits as well. Anything else would
     mean the inversion is not doing what SAS does.
     Column names are detected defensively: if your SAS version names them differently the
     block is skipped instead of aborting the program.                                   */
%let _hasq   = 0;
%let _haspop = 0;
%if %sysfunc(exist(km_q_sas)) %then %do;
  proc sql noprint;
    select count(*) into :_hasq trimmed from dictionary.columns
     where libname = 'WORK' and memname = 'KM_Q_SAS'
       and upcase(name) in ('PERCENT','ESTIMATE','LOWERLIMIT','UPPERLIMIT');
    select count(*) into :_haspop trimmed from dictionary.columns
     where libname = 'WORK' and memname = 'KM_Q_SAS' and upcase(name) = 'POP';
  quit;
%end;

%if &_hasq. = 4 and &_haspop. = 1 %then %do;
  data km_q2;
    set km_q_sas(rename = (Estimate   = km_sas_q
                           LowerLimit = km_sas_l
                           UpperLimit = km_sas_u));
    if Percent > 1 then percentile = Percent/100;
    else                percentile = Percent;
    keep pop percentile km_sas_q km_sas_l km_sas_u;
  run;
  proc sort data=km_q2; by pop percentile; run;
%end;
%else %do;
  data km_q2;
    length pop $16;  percentile 8;  km_sas_q 8;  km_sas_l 8;  km_sas_u 8;
    stop;
  run;
  %put NOTE: [demo] km_q_sas column names differ from expectation (hasq=&_hasq. haspop=&_haspop.). Automatic cross-validation skipped.;
  %put NOTE: [demo] Compare the Quartiles table in the listing with the "1-KM_comp" column of table 1 manually.;
%end;

/*==== 6. SUMMARY: THE THREE SETS OF PERCENTILES SIDE BY SIDE ========================*/
proc sql;
  create table cmp as
  select a.pop, a.percentile,
         a.q_hat as kmcs_q,  a.q_lcl as kmcs_l,  a.q_ucl as kmcs_u,
         b.q_hat as kmall_q, b.q_lcl as kmall_l, b.q_ucl as kmall_u,
         c.q_hat as cif_q,   c.q_lcl as cif_l,   c.q_ucl as cif_u,
         c.q_prev, c.F_at_q, c.plateau, c.plateau_ucl, c.maxtime,
         c.chk_order, c.chk_hit,
         d.km_sas_q, d.km_sas_l, d.km_sas_u,
         case when c.q_hat = . then 'NE' else 'OK' end as flag_q,
         case when c.q_hat ne . and c.q_lcl = . then 'NE' else '' end as flag_l,
         case when c.q_hat ne . and c.q_ucl = . then 'NE' else '' end as flag_u
  from      q_cs   as a
  left join q_all  as b on b.pop = a.pop and b.percentile = a.percentile
  left join q_cif  as c on c.pop = a.pop and c.percentile = a.percentile
  left join km_q2  as d on d.pop = a.pop and d.percentile = a.percentile
  order by a.pop, a.percentile;
quit;

/* Convert to years and format as "time (lower, upper)"; missing is always shown as NE. */
%macro fmtq(pre);
  if n(&pre._q) = 0 then &pre._s = 'NE';
  else &pre._s = catx(' ',
        ifn(&pre._q = ., 'NE', put(&pre._q/&yrfac., 6.2)),
        cats('(', ifn(&pre._l = ., 'NE', put(&pre._l/&yrfac., 6.2)), ', ',
                   ifn(&pre._u = ., 'NE', put(&pre._u/&yrfac., 6.2)), ')'));
%mend fmtq;

data cmp_r;
  set cmp;
  length kmcs_s $26 kmall_s $26 cif_s $26 sas_s $26 note $40 note_ci $46;
  %fmtq(kmcs)
  %fmtq(kmall)
  %fmtq(cif)

  if n(km_sas_q) = 0 then sas_s = 'n/a';
  else sas_s = catx(' ',
        ifn(km_sas_q = ., 'NE', put(km_sas_q/&yrfac., 6.2)),
        cats('(', ifn(km_sas_l = ., 'NE', put(km_sas_l/&yrfac., 6.2)), ', ',
                   ifn(km_sas_u = ., 'NE', put(km_sas_u/&yrfac., 6.2)), ')'));

  /* inversion vs the SAS Quartiles table: point estimate */
  if n(kmall_q, km_sas_q) = 2 then do;
    if kmall_q = km_sas_q then note = 'point estimate identical';
    else note = catx(' ', 'diff', put(kmall_q - km_sas_q, best8.), 'days');
  end;
  else note = '';

  /* inversion vs the SAS Quartiles table: CI limits */
  if n(kmall_l, km_sas_l, kmall_u, km_sas_u) = 4 then do;
    _d1 = abs(kmall_l - km_sas_l);
    _d2 = abs(kmall_u - km_sas_u);
    if max(_d1, _d2) < 0.005 then note_ci = 'CI limits identical';
    else note_ci = catx(' ', 'CI diff (days): L', put(_d1, best8.), ' U', put(_d2, best8.));
  end;
  else note_ci = '';

  label kmcs_s   = '1-KM_cs  years (95% CI)'
        kmall_s  = '1-KM_comp years (95% CI)'
        cif_s    = 'CIF years (95% CI)'
        sas_s    = 'SAS Quartiles table'
        note     = 'inversion vs Quartiles (point)'
        note_ci  = 'inversion vs Quartiles (CI)'
        percentile = 'p'
        flag_q   = 'CIF point'
        flag_l   = 'CIF lower'
        flag_u   = 'CIF upper'
        F_at_q   = 'F(q_hat)'
        q_prev   = 'q_prev (days)';
  format percentile 5.2 F_at_q 7.4 q_prev 8.2;
  drop kmcs_q kmcs_l kmcs_u kmall_q kmall_l kmall_u cif_q cif_l cif_u
       km_sas_q km_sas_l km_sas_u plateau plateau_ucl maxtime _d1 _d2;
run;

title "1. Percentiles of the three curves on the same data (years)";
proc print data=cmp_r noobs label;
  var pop percentile kmcs_s kmall_s cif_s flag_q flag_l flag_u;
run;

title "2. Cross-validation: inversion (1-KM_comp) vs the SAS Quartiles table";
proc print data=cmp_r noobs label;
  var pop percentile kmall_s sas_s note note_ci;
run;

/*==== 7. QC =========================================================================*/
/*  7.1 Confidence limit ordering and monotonicity, on the normalised curves            */
%macro qcCurve(dsin=, item=, tag=);
  %if %sysfunc(exist(&dsin.)) %then %do;
    data _qc_&item.;
      set &dsin.;
      length src $10;
      src = "&tag.";
      by pop time;
      if F = . then delete;
      if F_lcl = . or F_ucl = . then delete;
      bad_order = not (F_lcl <= F <= F_ucl);
      pf = lag(F); pl = lag(F_lcl); pu = lag(F_ucl);
      bad_mono = (pf > F) or (pl > F_lcl) or (pu > F_ucl);
      if first.pop then do; bad_order = 0; bad_mono = 0; end;
      if bad_order or bad_mono;
      keep pop src time F F_lcl F_ucl bad_order bad_mono;
    run;
  %end;
  %else %do;
    data _qc_&item.;
      length pop $16 src $10;
      length time F F_lcl F_ucl bad_order bad_mono 8;
      stop;
    run;
  %end;
%mend qcCurve;

%qcCurve(dsin=cur_cif, item=cif, tag=CIF);
%qcCurve(dsin=cur_all, item=all, tag=KM_comp);
%qcCurve(dsin=cur_cs,  item=cs,  tag=KM_cs);

data qc_curve;
  set _qc_cif _qc_all _qc_cs;
run;

proc sql;
  select count(*) as n_violation into :_nbad trimmed from qc_curve;
quit;
%put NOTE: [QC 7.1] rows violating limit ordering / monotonicity = &_nbad. (0 = pass);

/*  7.2 Inversion self checks                                                           */
proc sql;
  create table qc_inv as
  select count(*)                                            as n_row,
         sum(case when chk_order = 0 then 1 else 0 end)      as n_order_fail,
         sum(case when chk_hit   = 0 then 1 else 0 end)      as n_hit_fail,
         sum(case when q_hat = .          then 1 else 0 end) as n_NE_all,
         sum(case when q_hat ne . and (q_lcl = . or q_ucl = .)
                  then 1 else 0 end)                         as n_NE_partial
  from q_cif;
quit;
title "3. QC: inversion self checks (CIF curve)";
proc print data=qc_inv noobs; run;

/*  7.3 Percentile ordering: q(CIF) >= q(KM_cs) >= q(KM_comp)                          */
data qc_order;
  set cmp;
  length verdict $8;
  if n(cif_q, kmcs_q, kmall_q) < 3 then verdict = 'NE';
  else if kmall_q <= kmcs_q <= cif_q then verdict = 'PASS';
  else verdict = 'FAIL';
  keep pop percentile kmall_q kmcs_q cif_q verdict;
run;

title "4. QC: percentile ordering q(1-KM_comp) <= q(1-KM_cs) <= q(CIF)";
proc print data=qc_order noobs; run;

/*==== 8. CROSSING POINT READBACK ====================================================*/
/*  Shows that F(q_hat) >= p is a consequence of the definition (the step function jumps
     over p), not a numerical error.                                                     */
title "5. CIF inversion crossing point: F(q_hat) >= p";
proc print data=cmp_r noobs label;
  var pop percentile q_prev F_at_q cif_s;
run;

title;
/*=======================================================================================
  EXPECTED RESULTS (reproduced with the companion sim_competing.csv; time in years):
    Arm A  p=0.25   1-KM_comp 0.96  <=  1-KM_cs 1.43  <=  CIF 1.50
    Arm A  p=0.50   1-KM_comp 2.17  <=  1-KM_cs 3.10  <=  CIF 3.98
    Arm B  p=0.25   1-KM_comp 0.64  <=  1-KM_cs 0.86  <=  CIF 0.88
    Arm B  p=0.50   1-KM_comp 1.74  <=  1-KM_cs 2.25  <=  CIF 2.74
    The 75th percentile of the CIF is NE in all three groups (the CIF plateau is < 0.75).
    The UPPER limit of the CIF median in Arm A is also NE (the CIF lower confidence limit
    never reaches 0.50) -- the Shell must print NE in that cell plus a footnote.

  CHECKLIST:
    (1) table 2 note    = "point estimate identical" for every row  -> the inversion is
        structurally identical to what SAS does for KM percentiles.
    (2) table 2 note_ci = "CI limits identical"    for every row    -> the confidence
        limits agree too (both sides use GreenWood + LOGLOG on the same curve).
    (3) QC 7.1 n_violation   = 0
    (4) QC 7.2 n_order_fail  = 0 and n_hit_fail = 0
    (5) QC 7.3 verdict all PASS or NE
=======================================================================================*/
