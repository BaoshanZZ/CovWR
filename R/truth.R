estimate_truth_for_scenario <- function(scenario,
                                        n_pairs,
                                        chunk_size = 100000L,
                                        seed = NULL,
                                        params = default_dgm_params()) {
  if (!is.null(seed)) {
    set.seed(seed)
  }

  if (identical(as.character(scenario$effect), "null")) {
    return(data.table::data.table(
      scenario_id = scenario$scenario_id,
      truth_tau_w = NA_real_,
      truth_tau_l = NA_real_,
      truth_NB = 0,
      truth_LWR = 0,
      truth_LWO = 0
    ))
  }

  remaining <- as.integer(n_pairs)
  sum_w <- 0
  sum_l <- 0
  done <- 0L

  while (remaining > 0L) {
    n_chunk <- min(as.integer(chunk_size), remaining)
    treated <- generate_subject_data(
      N = n_chunk,
      A = 1L,
      effect = scenario$effect,
      prognosis = scenario$prognosis,
      outcome_shape = scenario$outcome_shape,
      effect_heterogeneity = scenario$effect_heterogeneity,
      endpoint_prognosis_pattern = scenario$endpoint_prognosis_pattern,
      params = params
    )
    control <- generate_subject_data(
      N = n_chunk,
      A = 0L,
      effect = scenario$effect,
      prognosis = scenario$prognosis,
      outcome_shape = scenario$outcome_shape,
      effect_heterogeneity = scenario$effect_heterogeneity,
      endpoint_prognosis_pattern = scenario$endpoint_prognosis_pattern,
      params = params
    )
    pairs <- hierarchical_pair_vectors(treated, control, params)
    sum_w <- sum_w + sum(pairs$win)
    sum_l <- sum_l + sum(pairs$loss)
    done <- done + n_chunk
    remaining <- remaining - n_chunk
  }

  tau_w <- sum_w / done
  tau_l <- sum_l / done
  measures <- calc_measures(tau_w, tau_l, params)

  data.table::data.table(
    scenario_id = scenario$scenario_id,
    truth_tau_w = tau_w,
    truth_tau_l = tau_l,
    truth_NB = unname(measures[["NB"]]),
    truth_LWR = unname(measures[["LWR"]]),
    truth_LWO = unname(measures[["LWO"]])
  )
}

estimate_truth_grid <- function(grid,
                                n_pairs,
                                chunk_size = 100000L,
                                seed = 20260617,
                                params = default_dgm_params()) {
  grid <- normalize_scenario_grid(grid)
  rows <- lapply(seq_len(nrow(grid)), function(i) {
    estimate_truth_for_scenario(
      scenario = grid[i],
      n_pairs = n_pairs,
      chunk_size = chunk_size,
      seed = seed + 1000L * i,
      params = params
    )
  })
  data.table::rbindlist(rows, use.names = TRUE)
}

truth_vector_for_scenario <- function(truth_row) {
  c(
    NB = truth_row$truth_NB,
    LWR = truth_row$truth_LWR,
    LWO = truth_row$truth_LWO
  )
}
