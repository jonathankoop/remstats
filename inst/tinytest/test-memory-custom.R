# memory = "custom": an arbitrary (tabulated) decay function of the lag
# --------------------------------------------------------------------
# The custom kernel is applied to the same statistics as memory = "decay"
# (everything derived from the weighted counts of past events), with the same
# reference time (the previous time point). So an exponential kernel tabulated
# on the observed lags has to reproduce memory = "decay" exactly.

# tie-oriented model, small edgelist (as in test-memory.R)
# --------------------------------------------------------------------
edgelist <- data.frame(
  time = 1:10,
  actor1 = c(1, 2, 1, 2, 3, 4, 2, 2, 2, 4),
  actor2 = c(3, 1, 3, 3, 2, 3, 1, 3, 4, 1)
)

reh <- remify(edgelist, model = "tie", riskset = "active")
effects <- ~ inertia() + reciprocity() +
  indegreeSender() + outdegreeReceiver() + totaldegreeDyad() +
  otp() + isp() + osp(scaling = "std") +
  psABBA() + recencyContinue() + rrankSend()

half_life <- 5
exp_kernel <- function(lag) exp(-lag * log(2) / half_life)
# all lags in this history are integers in 0:9
exp_table <- data.frame(lag = 0:10, weight = exp_kernel(0:10))

stats_decay <- remstats(
  reh = reh, tie_effects = effects,
  memory = "decay", memory_value = half_life, first = 1
)

# table as a data.frame
stats_custom <- remstats(
  reh = reh, tie_effects = effects,
  memory = "custom", memory_value = exp_table, first = 1
)
expect_equal(dimnames(stats_custom), dimnames(stats_decay))
for (sl in dimnames(stats_decay)[[3]]) {
  expect_equal(stats_custom[, , sl], stats_decay[, , sl], info = sl)
}

# table as a list
stats_list <- remstats(
  reh = reh, tie_effects = effects,
  memory = "custom", memory_value = list(lag = 0:10, weight = exp_kernel(0:10)),
  first = 1
)
expect_equal(unclass(stats_list), unclass(stats_custom))

# table as a matrix (unnamed columns: lag first, weight second; unsorted rows)
stats_matrix <- remstats(
  reh = reh, tie_effects = effects,
  memory = "custom", memory_value = cbind(10:0, exp_kernel(10:0)), first = 1
)
expect_equal(unclass(stats_matrix), unclass(stats_custom))

# a function of the lag (tabulated on a fine grid internally)
stats_fun <- remstats(
  reh = reh, tie_effects = effects,
  memory = "custom", memory_value = exp_kernel, first = 1
)
for (sl in dimnames(stats_decay)[[3]]) {
  expect_equal(stats_fun[, , sl], stats_decay[, , sl], tolerance = 1e-10, info = sl)
}

# subset of the history (reference time for the first row is the last time
# point before 'first')
stats_decay_sub <- remstats(
  reh = reh, tie_effects = effects,
  memory = "decay", memory_value = half_life, first = 4, last = 8
)
stats_custom_sub <- remstats(
  reh = reh, tie_effects = effects,
  memory = "custom", memory_value = exp_table, first = 4, last = 8
)
for (sl in dimnames(stats_decay_sub)[[3]]) {
  expect_equal(stats_custom_sub[, , sl], stats_decay_sub[, , sl], info = sl)
}

# a non-monotone step kernel, checked against a hand-computed inertia
# --------------------------------------------------------------------
step_kernel <- function(lag) ifelse(lag < 2, 1, ifelse(lag < 5, 0.5, 0))
step_table <- data.frame(lag = 0:10, weight = step_kernel(0:10))

stats_step <- remstats(
  reh = reh, tie_effects = ~ inertia(),
  memory = "custom", memory_value = step_table, first = 1
)
riskset <- attr(stats_step, "riskset")

# f(previous time point, event time): lag measured from the previous time point
f <- function(time, time_event) step_kernel(time - time_event)

inertia_step <- rbind(
  matrix(0, ncol = nrow(riskset)),
  c(f(1, 1), 0, 0, 0, 0, 0, 0),
  c(f(2, 1), f(2, 2), 0, 0, 0, 0, 0),
  c(f(3, 1) + f(3, 3), f(3, 2), 0, 0, 0, 0, 0),
  c(f(4, 1) + f(4, 3), f(4, 2), f(4, 4), 0, 0, 0, 0),
  c(f(5, 1) + f(5, 3), f(5, 2), f(5, 4), 0, f(5, 5), 0, 0),
  c(f(6, 1) + f(6, 3), f(6, 2), f(6, 4), 0, f(6, 5), 0, f(6, 6)),
  c(f(7, 1) + f(7, 3), f(7, 2) + f(7, 7), f(7, 4), 0, f(7, 5), 0, f(7, 6)),
  c(f(8, 1) + f(8, 3), f(8, 2) + f(8, 7), f(8, 4) + f(8, 8), 0, f(8, 5), 0, f(8, 6)),
  c(f(9, 1) + f(9, 3), f(9, 2) + f(9, 7), f(9, 4) + f(9, 8), f(9, 9), f(9, 5), 0, f(9, 6))
)
expect_equal(stats_step[, , "inertia"], inertia_step)

# nearest-neighbour lookup: lags between grid points take the closest grid
# weight (ties go to the smaller lag), lags beyond the grid the last weight
coarse_table <- data.frame(lag = c(0, 3), weight = c(1, 0.25))
stats_coarse <- remstats(
  reh = reh, tie_effects = ~ inertia(),
  memory = "custom", memory_value = coarse_table, first = 1
)
# lag 0, 1 -> 1; lag 2, 3, ... -> 0.25 (lag 1.5 would tie and go to lag 0)
nn_kernel <- function(lag) ifelse(lag <= 1, 1, 0.25)
f <- function(time, time_event) nn_kernel(time - time_event)
inertia_nn <- rbind(
  matrix(0, ncol = nrow(riskset)),
  c(f(1, 1), 0, 0, 0, 0, 0, 0),
  c(f(2, 1), f(2, 2), 0, 0, 0, 0, 0),
  c(f(3, 1) + f(3, 3), f(3, 2), 0, 0, 0, 0, 0),
  c(f(4, 1) + f(4, 3), f(4, 2), f(4, 4), 0, 0, 0, 0),
  c(f(5, 1) + f(5, 3), f(5, 2), f(5, 4), 0, f(5, 5), 0, 0),
  c(f(6, 1) + f(6, 3), f(6, 2), f(6, 4), 0, f(6, 5), 0, f(6, 6)),
  c(f(7, 1) + f(7, 3), f(7, 2) + f(7, 7), f(7, 4), 0, f(7, 5), 0, f(7, 6)),
  c(f(8, 1) + f(8, 3), f(8, 2) + f(8, 7), f(8, 4) + f(8, 8), 0, f(8, 5), 0, f(8, 6)),
  c(f(9, 1) + f(9, 3), f(9, 2) + f(9, 7), f(9, 4) + f(9, 8), f(9, 9), f(9, 5), 0, f(9, 6))
)
expect_equal(stats_coarse[, , "inertia"], inertia_nn)

# undirected events
# --------------------------------------------------------------------
reh_und <- remify(edgelist, model = "tie", directed = FALSE)
effects_und <- ~ inertia() + degreeMin() + degreeMax() + totaldegreeDyad() + sp()

stats_und_decay <- remstats(
  reh = reh_und, tie_effects = effects_und,
  memory = "decay", memory_value = half_life, first = 1
)
stats_und_custom <- remstats(
  reh = reh_und, tie_effects = effects_und,
  memory = "custom", memory_value = exp_table, first = 1
)
for (sl in dimnames(stats_und_decay)[[3]]) {
  expect_equal(stats_und_custom[, , sl], stats_und_decay[, , sl], info = sl)
}

# case-control sampling: sampled custom == sampled decay (exponential table)
# --------------------------------------------------------------------
data(history, package = "remstats")
data(info, package = "remstats")
reh_hist <- remify(edgelist = history, model = "tie", riskset = "active")

samp_effects <- ~ inertia() + reciprocity() +
  indegreeSender(scaling = "prop") + outdegreeReceiver() + totaldegreeDyad() +
  otp() + isp() + osp(scaling = "std") + psABBA() + recencyContinue()

lags_hist <- with(reh_hist$edgelist, as.numeric(outer(time, time, "-")))
lags_hist <- sort(unique(lags_hist[lags_hist >= 0]))
half_life_hist <- 1000
exp_table_hist <- data.frame(
  lag = lags_hist, weight = exp(-lags_hist * log(2) / half_life_hist)
)

samp_decay <- tomstats(
  samp_effects, reh = reh_hist, memory = "decay", memory_value = half_life_hist,
  sampling = TRUE, samp_num = 5L, seed = 1L, first = 2, last = 30
)
samp_custom <- tomstats(
  samp_effects, reh = reh_hist, memory = "custom", memory_value = exp_table_hist,
  sampling = TRUE, samp_num = 5L, seed = 1L, first = 2, last = 30
)
expect_equal(attr(samp_custom, "sample_map"), attr(samp_decay, "sample_map"))
for (sl in dimnames(samp_decay)[[3]]) {
  expect_equal(samp_custom[, , sl], samp_decay[, , sl], tolerance = 1e-12, info = sl)
}

# case-control sampling: sampled == full at the sampled dyads for a smooth
# non-exponential kernel (all inertia-derived statistics)
# --------------------------------------------------------------------
max_lag_hist <- diff(range(reh_hist$edgelist$time))
smooth_kernel <- function(lag) {
  (1 + cos(pi * pmin(lag, max_lag_hist) / max_lag_hist)) / 2 * (1 + sin(lag / 500)) / 2
}
smooth_table <- data.frame(
  lag = seq(0, max_lag_hist, length.out = 10000),
  weight = smooth_kernel(seq(0, max_lag_hist, length.out = 10000))
)

check_sampled_equals_full <- function(effects, memory_value, samp_num = 5L,
                                      seed = 1L, tol = 1e-12) {
  ts_samp <- tomstats(
    effects, reh = reh_hist, memory = "custom", memory_value = memory_value,
    sampling = TRUE, samp_num = samp_num, seed = seed, first = 2, last = 40
  )
  ts_full <- tomstats(
    effects, reh = reh_hist, memory = "custom", memory_value = memory_value,
    sampling = FALSE, first = 2, last = 40
  )
  sample_map <- attr(ts_samp, "sample_map")
  expect_true(!is.null(sample_map))
  for (m in seq_len(dim(ts_samp)[1])) {
    for (s in seq_len(dim(ts_samp)[2])) {
      expect_equal(
        as.numeric(ts_samp[m, s, ]),
        as.numeric(ts_full[m, sample_map[m, s], ]),
        tolerance = tol,
        info = sprintf("m=%d s=%d", m, s)
      )
    }
  }
  invisible(TRUE)
}

check_sampled_equals_full(~ inertia() + reciprocity(), smooth_table)
check_sampled_equals_full(
  ~ inertia(scaling = "prop") + reciprocity(scaling = "prop"), smooth_table
)
check_sampled_equals_full(
  ~ indegreeSender() + outdegreeSender() + indegreeReceiver() + outdegreeReceiver() +
    totaldegreeSender() + totaldegreeReceiver(), smooth_table
)
check_sampled_equals_full(
  ~ indegreeSender(scaling = "prop") + outdegreeReceiver(scaling = "prop"),
  smooth_table
)
check_sampled_equals_full(~ totaldegreeDyad(), smooth_table)
check_sampled_equals_full(~ otp() + itp() + isp() + osp(), smooth_table)
# (scaling = "std" standardizes over the sampled dyads only, so sampled and
# full differ under every memory type; not tested here, as in
# test_compare_sampled_stats4.R)
# the same with the kernel passed as a function
check_sampled_equals_full(~ inertia() + otp(), smooth_kernel)

# the smooth kernel is not exponential: custom must differ from decay
stats_smooth <- tomstats(
  ~ inertia(), reh = reh_hist, memory = "custom", memory_value = smooth_table,
  first = 2, last = 40
)
stats_exp <- tomstats(
  ~ inertia(), reh = reh_hist, memory = "decay", memory_value = half_life_hist,
  first = 2, last = 40
)
expect_false(isTRUE(all.equal(stats_smooth[, , "inertia"], stats_exp[, , "inertia"])))

# actor-oriented model: custom (exponential table) == decay
# --------------------------------------------------------------------
reh_actor <- remify(edgelist, model = "actor")
sender_effects <- ~ indegreeSender() + outdegreeSender() + totaldegreeSender() +
  recencySendSender() + psABA()
receiver_effects <- ~ inertia() + reciprocity() +
  indegreeReceiver() + outdegreeReceiver() + totaldegreeReceiver() +
  otp() + itp() + osp() + isp() +
  rrankSend() + recencyContinue() + psABBA()

aom_decay <- remstats(
  reh = reh_actor, sender_effects = sender_effects,
  receiver_effects = receiver_effects,
  memory = "decay", memory_value = half_life, first = 1
)
aom_custom <- remstats(
  reh = reh_actor, sender_effects = sender_effects,
  receiver_effects = receiver_effects,
  memory = "custom", memory_value = exp_table, first = 1
)
expect_equal(unclass(aom_custom$sender_stats), unclass(aom_decay$sender_stats))
expect_equal(unclass(aom_custom$receiver_stats), unclass(aom_decay$receiver_stats))

aom_fun <- remstats(
  reh = reh_actor, sender_effects = sender_effects,
  receiver_effects = receiver_effects,
  memory = "custom", memory_value = exp_kernel, first = 1
)
expect_equal(unclass(aom_fun$sender_stats), unclass(aom_decay$sender_stats),
  tolerance = 1e-10)
expect_equal(unclass(aom_fun$receiver_stats), unclass(aom_decay$receiver_stats),
  tolerance = 1e-10)

# proportional scaling (uses the weighted degree of the sender)
aom_prop_decay <- remstats(
  reh = reh_actor, receiver_effects = ~ inertia(scaling = "prop"),
  memory = "decay", memory_value = half_life, first = 1
)
aom_prop_custom <- remstats(
  reh = reh_actor, receiver_effects = ~ inertia(scaling = "prop"),
  memory = "custom", memory_value = exp_table, first = 1
)
expect_equal(unclass(aom_prop_custom$receiver_stats),
  unclass(aom_prop_decay$receiver_stats))

# actor-oriented step kernel: receiver inertia equals the tie-oriented inertia
# of the sender's dyads
aom_step <- remstats(
  reh = reh_actor, receiver_effects = ~ inertia(),
  memory = "custom", memory_value = step_table, first = 1
)
tie_step <- remstats(
  reh = remify(edgelist, model = "tie"), tie_effects = ~ inertia(),
  memory = "custom", memory_value = step_table, first = 1
)
tie_riskset <- attr(tie_step, "riskset")
actor_names <- as.character(reh_actor$meta$dictionary$actors[, 1])
for (m in seq_len(nrow(edgelist))) {
  s <- edgelist$actor1[m]
  cols <- which(as.character(tie_riskset$sender) == as.character(s))
  j <- match(as.character(tie_riskset$receiver[cols]), actor_names)
  expect_equal(
    as.numeric(aom_step$receiver_stats[m, j, "inertia"]),
    as.numeric(tie_step[m, cols, "inertia"]),
    info = paste("row", m)
  )
}

# invalid input
# --------------------------------------------------------------------
expect_error(
  remstats(reh = reh, tie_effects = ~ inertia(), memory = "custom"),
  "function of the lag or a table"
)
expect_error(
  remstats(reh = reh, tie_effects = ~ inertia(), memory = "custom",
    memory_value = 5),
  "function of the lag or a table"
)
expect_error(
  remstats(reh = reh, tie_effects = ~ inertia(), memory = "custom",
    memory_value = data.frame(lag = c(-1, 0, 1), weight = c(1, 1, 1))),
  "non-negative"
)
expect_error(
  remstats(reh = reh, tie_effects = ~ inertia(), memory = "custom",
    memory_value = data.frame(lag = c(0, 0, 1), weight = c(1, 1, 1))),
  "unique"
)
expect_error(
  remstats(reh = reh, tie_effects = ~ inertia(), memory = "custom",
    memory_value = data.frame(lag = c(0, 1), weight = c(1, NA))),
  "finite weight"
)
expect_error(
  remstats(reh = reh, tie_effects = ~ inertia(), memory = "custom",
    memory_value = data.frame(lags = 0:1, w = c(1, 1))),
  "'lag' and 'weight'"
)
expect_error(
  remstats(reh = reh, tie_effects = ~ inertia(), memory = "custom",
    memory_value = cbind(0:1, c(1, 1), c(1, 1))),
  "two columns"
)
