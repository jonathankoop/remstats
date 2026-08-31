#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# dev-bench-sampled.R
#
# Where does the case-control sampled path actually fall over?
#
#   Rscript dev-bench-sampled.R              # full sweep (~minutes)
#   Rscript dev-bench-sampled.R --quick      # small grid, ~30s
#   Rscript dev-bench-sampled.R --profile    # + Rprof on one config
#
# The hypothesis under test: sampling is supposed to make cost independent of
# the risk set, but tomstats(sampling = TRUE) still calls prepare_tomstats(),
# which materialises
#
#     prepR         : D x 4 double matrix,        D = N(N-1)C
#     risksetMatrix : N x (N*C) double matrix
#
# both O(N^2). If that is the wall, it is an R-side / remify problem and no
# amount of work in tomstats_sampled.cpp helps. This script separates:
#
#   t_remify   remify()                       -- edgelist prep
#   t_prepare  prepare_tomstats()             -- the O(N^2) construction
#   t_rest     tomstats(sampling=TRUE) minus the above
#              = R-side sampling machinery + compute_stats_tie_sampled()
#
# and reports peak memory + the size of the two N^2 objects alongside.
#
# Results stream to bench-sampled-results.csv as they are produced, so a
# config that dies by allocation failure still leaves everything before it.
# ---------------------------------------------------------------------------

suppressMessages({
  library(remify)
  library(remstats)
})

args    <- commandArgs(trailingOnly = TRUE)
QUICK   <- "--quick"   %in% args
PROFILE <- "--profile" %in% args
BIG     <- "--big"     %in% args
OUT     <- "bench-sampled-results.csv"

# Refuse to build a risk set projected larger than this (Mb). The point of the
# sweep is to find the wall, not to drive the machine into swap discovering it.
# --big raises the cap to 16 GB if you want to watch it actually die.
CAP_MB <- if (BIG) 16 * 1024 else 3 * 1024

# Abandon a sweep once a single config crosses this (seconds). Compiled code
# ignores setTimeLimit(), so this is checked between configs, not within one.
BUDGET <- if (QUICK) 20 else 120

# ---------------------------------------------------------------------------
# effect sets
#
# "baseline" isolates pure fixed overhead: no statistic is computed at all, so
# whatever it costs is setup. "core" adds the stats whose emit path is already
# sampling-proportional. "triad" adds the one stat with per-cell work in N.
# "wide" is core repeated to 12 effects -- the difference between core and wide
# is per-effect fixed cost (time_points sort, risksetMatrix.max(), D-sized
# lookup tables), which is the thing worth measuring.
# ---------------------------------------------------------------------------
EFFECT_SETS <- list(
  baseline = ~ baseline(),
  core     = ~ inertia() + reciprocity() + indegreeSender() + outdegreeSender(),
  wide     = ~ inertia() + reciprocity() +
               indegreeSender() + indegreeReceiver() +
               outdegreeSender() + outdegreeReceiver() +
               totaldegreeSender() + totaldegreeReceiver() +
               degreeMin() + degreeMax() + degreeDiff() + totaldegreeDyad(),
  triad    = ~ inertia() + reciprocity() + otp() + isp()
)

# ---------------------------------------------------------------------------
# synthetic data
#
# Preferential-attachment-ish sender/receiver draws: uniform actors would give
# a network with no inertia and no triangles, which flatters every sparse
# optimisation. alpha controls skew; alpha = 0 is uniform.
# ---------------------------------------------------------------------------
make_history <- function(N, E, alpha = 0.8, seed = 42) {
  set.seed(seed)
  w  <- (seq_len(N))^(-alpha)
  w  <- w / sum(w)
  a1 <- sample.int(N, E, replace = TRUE, prob = w)
  a2 <- sample.int(N, E, replace = TRUE, prob = w)
  # no self-loops
  bad <- which(a1 == a2)
  while (length(bad)) {
    a2[bad] <- sample.int(N, length(bad), replace = TRUE, prob = w)
    bad <- bad[a1[bad] == a2[bad]]
  }
  data.frame(
    time   = seq_len(E),                    # distinct times => M == E under "pt"
    actor1 = as.character(a1),
    actor2 = as.character(a2),
    stringsAsFactors = FALSE
  )
}

# ---------------------------------------------------------------------------
# measurement helpers
# ---------------------------------------------------------------------------
peak_mb <- function() sum(gc(verbose = FALSE)[, 6])

timed <- function(expr) {
  t0 <- proc.time()[["elapsed"]]
  val <- force(expr)
  list(value = val, secs = proc.time()[["elapsed"]] - t0)
}

# predicted size of the two O(N^2) objects, in Mb, without building them
predict_mb <- function(N, C = 1) {
  D  <- N * (N - 1) * C
  c(prepR_mb  = D * 4 * 8 / 2^20,          # D x 4 doubles
    rsmat_mb  = N * N * C * 8 / 2^20)      # N x (N*C) doubles
}

emit <- function(row) {
  write.table(as.data.frame(row), OUT,
              sep = ",", row.names = FALSE,
              col.names = !file.exists(OUT),
              append = file.exists(OUT))
  invisible(row)
}

# ---------------------------------------------------------------------------
# one configuration
# ---------------------------------------------------------------------------
run_one <- function(N, E, S, set_name, memory = "full", memory_value = NA) {

  eff  <- EFFECT_SETS[[set_name]]
  pred <- predict_mb(N)

  cat(sprintf("N=%-7d E=%-8d S=%-3d %-9s ", N, E, S, set_name))
  utils::flush.console()

  base <- list(
    N = N, E = E, S = S, effects = set_name, memory = memory,
    n_eff = length(attr(terms(eff), "term.labels")),
    pred_prepR_mb = round(pred[["prepR_mb"]], 1),
    pred_rsmat_mb = round(pred[["rsmat_mb"]], 1)
  )

  fail <- function(where, e) {
    cat(sprintf("FAILED (%s: %s)\n", where, conditionMessage(e)))
    emit(c(base, list(t_remify = NA, t_prepare = NA, t_total = NA, t_rest = NA,
                      prepR_mb = NA, rsmat_mb = NA, peak_mb = NA,
                      status = paste0("fail:", where))))
    NULL
  }

  if (sum(pred) > CAP_MB) {
    cat(sprintf("SKIPPED (projected %.1f GB of risk set; --big to force)\n",
                sum(pred) / 1024))
    emit(c(base, list(t_remify = NA, t_prepare = NA, t_total = NA, t_rest = NA,
                      prepR_mb = NA, rsmat_mb = NA, peak_mb = NA,
                      status = "skipped:projected_size")))
    return(invisible(NULL))
  }

  el <- make_history(N, E)
  gc(reset = TRUE, verbose = FALSE)

  reh <- tryCatch(
    timed(remify::remify(edgelist = el, model = "tie",
                         actors = as.character(seq_len(N)),
                         riskset = "full", directed = TRUE)),
    error = function(e) fail("remify", e))
  if (is.null(reh)) return(invisible(NULL))

  # the suspected wall, timed on its own
  prep <- tryCatch(
    timed(remstats:::prepare_tomstats(
      effects = eff, reh = reh$value,
      memory = memory, memory_value = memory_value,
      start = 2, stop = Inf, method = "pt")),
    error = function(e) fail("prepare_tomstats", e))
  if (is.null(prep)) return(invisible(NULL))

  prepR_mb  <- as.numeric(object.size(prep$value$riskset))       / 2^20
  rsmat_mb  <- as.numeric(object.size(prep$value$risksetMatrix)) / 2^20
  prep_secs <- prep$secs          # keep before dropping the inputs list
  rm(prep); gc(verbose = FALSE)

  tot <- tryCatch(
    timed(tomstats(eff, reh = reh$value,
                   memory = memory, memory_value = memory_value,
                   sampling = TRUE, samp_num = as.integer(S), seed = 1L,
                   first = 2, last = Inf)),
    error = function(e) fail("tomstats", e))
  if (is.null(tot)) return(invisible(NULL))

  t_rest <- tot$secs - prep_secs
  pk <- peak_mb()

  cat(sprintf("| remify %6.2fs  prepare %6.2fs  rest %6.2fs  | prepR %7.1fMb  rsmat %7.1fMb  peak %7.1fMb\n",
              reh$secs, prep_secs, t_rest, prepR_mb, rsmat_mb, pk))

  emit(c(base, list(
    t_remify = round(reh$secs, 3), t_prepare = round(prep_secs, 3),
    t_total  = round(tot$secs, 3), t_rest    = round(t_rest, 3),
    prepR_mb = round(prepR_mb, 1), rsmat_mb  = round(rsmat_mb, 1),
    peak_mb  = round(pk, 1), status = "ok")))

  invisible(list(secs = tot$secs))
}

# ---------------------------------------------------------------------------
# sweeps
# ---------------------------------------------------------------------------
if (file.exists(OUT)) file.remove(OUT)

Ns <- if (QUICK) {
  c(100, 250, 500)
} else if (BIG) {
  c(100, 250, 500, 1000, 2000, 4000, 8000, 16000)
} else {
  c(100, 250, 500, 1000, 2000, 4000)
}
Es <- if (QUICK) c(2000)          else c(2000, 10000, 50000)
S  <- 20L

cat("\n== sweep 1: N (actors) at fixed E, effect set 'core' ==\n")
cat("   isolates the O(N^2) setup; S and E are held constant, so anything\n")
cat("   that grows here is risk-set construction, not statistic computation\n\n")
for (N in Ns) {
  r <- run_one(N, E = Es[1], S = S, set_name = "core")
  if (!is.null(r) && r$secs > BUDGET) { cat("   (budget exceeded, stopping sweep)\n"); break }
}

cat("\n== sweep 2: E (events) at fixed N=500, effect set 'core' ==\n")
cat("   the part that should scale linearly if the emit path is healthy\n\n")
for (E in Es) {
  r <- run_one(500, E = E, S = S, set_name = "core")
  if (!is.null(r) && r$secs > BUDGET) { cat("   (budget exceeded, stopping sweep)\n"); break }
}

cat("\n== sweep 3: fixed vs per-effect cost at N=1000, E=10000 ==\n")
cat("   baseline computes nothing, so its t_rest is pure overhead.\n")
cat("   (wide - core) / 8 extra effects = marginal cost per effect;\n")
cat("   compare against core/4. If they are similar, per-effect fixed setup\n")
cat("   dominates and the shared-context refactor is the win.\n\n")
for (set_name in c("baseline", "core", "wide", "triad")) {
  r <- run_one(1000, E = 10000, S = S, set_name = set_name)
  if (!is.null(r) && r$secs > BUDGET) cat("   (over budget; larger sets skipped)\n")
}

cat("\n== sweep 4: S (sampled dyads) at N=1000, E=10000 ==\n")
cat("   cost should be ~flat in S if fixed overhead dominates,\n")
cat("   ~linear in S if the emit path does\n\n")
for (s in c(5L, 20L, 50L, 100L)) run_one(1000, E = 10000, S = s, set_name = "core")

cat("\n== sweep 5: memory kinds at N=1000, E=10000 ==\n\n")
run_one(1000, 10000, S, "core", memory = "full")
run_one(1000, 10000, S, "core", memory = "decay",  memory_value = 100)
run_one(1000, 10000, S, "core", memory = "window", memory_value = 100)

# ---------------------------------------------------------------------------
# projection past the point where it dies
# ---------------------------------------------------------------------------
cat("\n== projected size of the two O(N^2) objects (C = 1) ==\n\n")
fmt_mb <- function(mb) {
  if (mb < 1024) sprintf("%.1f MB", mb)
  else if (mb < 1024^2) sprintf("%.1f GB", mb / 1024)
  else sprintf("%.1f TB", mb / 1024^2)
}
proj <- do.call(rbind, lapply(c(1e3, 1e4, 1e5, 1e6), function(N) {
  p <- predict_mb(N)
  data.frame(N = format(N, scientific = FALSE, big.mark = ","),
             prepR = fmt_mb(p[["prepR_mb"]]),
             risksetMatrix = fmt_mb(p[["rsmat_mb"]]))
}))
print(proj, row.names = FALSE)

if (PROFILE) {
  cat("\n== Rprof: N=2000, E=10000, core ==\n\n")
  el   <- make_history(2000, 10000)
  reh  <- remify::remify(edgelist = el, model = "tie",
                         actors = as.character(seq_len(2000)),
                         riskset = "full", directed = TRUE)
  Rprof("bench-sampled.Rprof", interval = 0.005, memory.profiling = TRUE)
  invisible(tomstats(EFFECT_SETS$core, reh = reh, memory = "full",
                     sampling = TRUE, samp_num = 20L, seed = 1L,
                     first = 2, last = Inf))
  Rprof(NULL)
  print(head(summaryRprof("bench-sampled.Rprof")$by.self, 20))
  cat("\n(compiled time shows as the .Call frame; R-side hotspots are the\n")
  cat(" paste()/match() key building in tomstats.R and prepare_tomstats)\n")
}

cat(sprintf("\nWrote %s\n", OUT))
