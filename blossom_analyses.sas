/*========================================================
  Analysis plan
  C0 / T0: Primary - visualization, task type, interaction.
  C1 / T1: Sensitivity 1 - add categorical block order.
  C2 / T2: Sensitivity 2 - also add vis*vis_order.
  All models use a participant random intercept only.
  Visualization comparisons use Tukey adjustment within
  each task, separately for correctness and completion time.
  D1-D3 retain the block-adjusted transformation diagnostics.
  Estimation and DF settings are retained from this script;
  C0/T0 restore the original fixed/random model structure.
========================================================*/
title;
title "Study data - Import";
proc import file="blossom_correctness_time.csv"
    DBMS=CSV
	OUT=BLOSSOM.CORRECTNESS;
	GETNAMES=YES;
run;

/* Optional: check variable names and levels */
title "Study data - Variable information";
proc contents data=BLOSSOM.CORRECTNESS; run;
title "Study data - Visualization, task, and correctness frequencies";
proc freq data=BLOSSOM.CORRECTNESS;
    tables vis task_order correctness / missing;
run;

/* Correctness: primary model with original fixed-effects structure */
title "C0: Correctness - Primary model (visualization, task, and interaction)";
proc glimmix data=BLOSSOM.CORRECTNESS method=laplace;
    class participant_id
          vis (ref="blossom")
          task_order (ref="1");

    model correctness(event='1') =
          vis
          task_order
          vis*task_order
          / dist=binary link=logit solution ddfm=residual;

    random intercept / subject=participant_id;

    /* Estimated means only; no overall pairwise comparisons */
    lsmeans vis / ilink cl;
    lsmeans task_order / ilink cl;

    /* Three visualization comparisons within each task */
    lsmeans vis*task_order / ilink slicediff=task_order
                            oddsratio cl adjust=tukey;
run;

/* Correctness: sensitivity 1, adjusting for visualization block order */
title "C1: Correctness - Sensitivity 1 (add visualization block order)";
proc glimmix data=BLOSSOM.CORRECTNESS method=laplace;
    class participant_id
          vis (ref="blossom")
          task_order (ref="1")
          vis_order (ref="1");

    model correctness(event='1') =
          vis
          task_order
          vis*task_order
          vis_order
          / dist=binary link=logit solution ddfm=residual;

    random intercept / subject=participant_id;

    /* Estimated means only; no overall pairwise comparisons */
    lsmeans vis / ilink cl;
    lsmeans task_order / ilink cl;
    lsmeans vis*task_order / ilink slicediff=task_order oddsratio cl adjust=tukey;
run;

/* Correctness: sensitivity 2, allowing visualization-specific block effects */
title "C2: Correctness - Sensitivity 2 (add visualization by block order interaction)";
proc glimmix data=BLOSSOM.CORRECTNESS method=laplace;
    class participant_id
          vis (ref="blossom")
          task_order (ref="1")
          vis_order (ref="1");

    model correctness(event='1') =
          vis
          task_order
          vis*task_order
          vis_order
          vis*vis_order
          / dist=binary link=logit solution ddfm=residual;

    random intercept / subject=participant_id;

    /* Estimated correctness for each visualization in each block */
    lsmeans vis*vis_order / ilink cl;

    /* Three visualization comparisons within each task */
    lsmeans vis*task_order / ilink slicediff=task_order
                            oddsratio cl adjust=tukey;
run;


/*========================================================
  1. Create transformed variables from imported dataset
========================================================*/
data BLOSSOM.TIME_DATA;
    set BLOSSOM.CORRECTNESS;

    /* Raw time */
    raw_duration = duration_ms;

    /* Log transform */
    if duration_ms > 0 then log_duration = log(duration_ms);
    else log_duration = .;

    /* Box-Cox transform with lambda = -0.341 */
    if duration_ms > 0 then bc_duration = (duration_ms**(-0.341) - 1)/(-0.341);
    else bc_duration = .;
run;

/* Optional check */
title "Completion time - Descriptive checks of transformed variables";
proc means data=BLOSSOM.TIME_DATA n nmiss min max;
    var raw_duration log_duration bc_duration;
run;

/*========================================================
  2. D1: Raw duration diagnostics (block-adjusted model)
========================================================*/
ods output FitStatistics=Fit_Raw
           Tests3=Tests_Raw
           SliceDiffs=SliceDiffs_Raw
           LSMeans=LSMeans_Raw;

title "D1: Completion time diagnostics - Raw duration (block-adjusted model)";
proc mixed data=BLOSSOM.TIME_DATA method=reml plots=studentpanel;
    class participant_id vis task_order vis_order;
    model raw_duration = vis task_order vis*task_order vis_order / solution ddfm=kr;
    random intercept / subject=participant_id;
    lsmeans vis*task_order / cl;
    slice vis*task_order / sliceby=task_order diff
                          cl adjust=tukey adjdfe=row;
run;

/*========================================================
  3. D2: Log duration diagnostics (block-adjusted model)
========================================================*/
ods output FitStatistics=Fit_Log
           Tests3=Tests_Log
           SliceDiffs=SliceDiffs_Log
           LSMeans=LSMeans_Log;

title "D2: Completion time diagnostics - Log duration (block-adjusted model)";
proc mixed data=BLOSSOM.TIME_DATA method=reml plots=studentpanel;
    class participant_id vis task_order vis_order;
    model log_duration = vis task_order vis*task_order vis_order / solution ddfm=kr;
    random intercept / subject=participant_id;
    lsmeans vis*task_order / cl;
    slice vis*task_order / sliceby=task_order diff
                          cl adjust=tukey adjdfe=row;
run;

/*========================================================
  4. D3: Box-Cox duration diagnostics (block-adjusted model)
========================================================*/
ods output FitStatistics=Fit_BoxCox
           Tests3=Tests_BoxCox
           SliceDiffs=SliceDiffs_BoxCox
           LSMeans=LSMeans_BoxCox;

title "D3: Completion time diagnostics - Box-Cox duration (block-adjusted model)";
proc mixed data=BLOSSOM.TIME_DATA method=reml plots=studentpanel;
    class participant_id vis task_order vis_order;
    model bc_duration = vis task_order vis*task_order vis_order / solution ddfm=kr;
    random intercept / subject=participant_id;
    lsmeans vis*task_order / cl;
    slice vis*task_order / sliceby=task_order diff
                          cl adjust=tukey adjdfe=row;
run;

/*========================================================
  5. Label the fit statistics
========================================================*/
data Fit_Raw;
    set Fit_Raw;
    Model = "Raw Duration";
run;

data Fit_Log;
    set Fit_Log;
    Model = "Log Duration";
run;

data Fit_BoxCox;
    set Fit_BoxCox;
    Model = "BoxCox Duration";
run;

/*========================================================
  6. Combine fit statistics for comparison
========================================================*/
data Fit_Compare;
    set Fit_Raw Fit_Log Fit_BoxCox;
run;

title "D4: Completion time diagnostics - Fit statistics by transformation";
proc print data=Fit_Compare;
run;

/*========================================================
  7. Label and combine Type III tests
========================================================*/
data Tests_Raw;
    set Tests_Raw;
    Model = "Raw Duration";
run;

data Tests_Log;
    set Tests_Log;
    Model = "Log Duration";
run;

data Tests_BoxCox;
    set Tests_BoxCox;
    Model = "BoxCox Duration";
run;

data Tests_Compare;
    set Tests_Raw Tests_Log Tests_BoxCox;
run;

title "D5: Completion time diagnostics - Fixed effects by transformation";
proc print data=Tests_Compare;
run;

/* Completion time: primary model with original fixed-effects structure */
title "T0: Completion time - Primary model (Box-Cox transformed)";
proc glimmix data=BLOSSOM.TIME_DATA method=rspl;
    class participant_id
          vis (ref="blossom")
          task_order (ref="1");

    model bc_duration =
          vis
          task_order
          vis*task_order
          / dist=normal link=identity solution ddfm=kr;

    random intercept / subject=participant_id;

    lsmeans vis*task_order / slicediff=task_order
                            cl adjust=tukey adjdfe=row;
run;

/* Completion time: sensitivity 1, adjusting for visualization block order */
title "T1: Completion time - Sensitivity 1 (add visualization block order)";
proc glimmix data=BLOSSOM.TIME_DATA method=rspl;
    class participant_id
          vis (ref="blossom")
          task_order (ref="1")
          vis_order (ref="1");

    model bc_duration =
          vis
          task_order
          vis*task_order
          vis_order
          / dist=normal link=identity solution ddfm=kr;

    random intercept / subject=participant_id;

    lsmeans vis*task_order / slicediff=task_order
                            cl adjust=tukey adjdfe=row;
run;

/* Completion time: sensitivity 2, allowing visualization-specific block effects */
title "T2: Completion time - Sensitivity 2 (add visualization by block order interaction)";
proc glimmix data=BLOSSOM.TIME_DATA method=rspl;
    class participant_id
          vis (ref="blossom")
          task_order (ref="1")
          vis_order (ref="1");

    model bc_duration =
          vis
          task_order
          vis*task_order
          vis_order
          vis*vis_order
          / dist=normal link=identity solution ddfm=kr;

    random intercept / subject=participant_id;

    /* Estimated means for each visualization in each block */
    lsmeans vis*vis_order / cl;

    /* Visualization comparisons within each task */
    lsmeans vis*task_order / slicediff=task_order
                            cl adjust=tukey adjdfe=row;
run;

/* Clear the final model title before any later SAS output. */
title;
