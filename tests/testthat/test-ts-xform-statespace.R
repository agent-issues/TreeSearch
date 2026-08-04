# Tier 2: skipped on CRAN; see tests/testing-strategy.md
skip_on_cran()

# Regressions for red-team T-393, T-394 and T-401: how a hierarchy block's
# secondary state space is derived, what happens when a secondary has no
# observed state, and what `MaximizeParsimony()` may score its own output with.

library("TreeTools")

# `MatrixToPhyDat()` rather than the `phangorn::phyDat()` helper the rest of the
# xform suite uses: only it puts an ambiguity token such as "{01}" in the
# contrast matrix as the state SET it denotes, which is the whole subject here.
xss_dat <- function(mat) MatrixToPhyDat(mat)

xss_tree <- function() ape::read.tree(text = "(((t1,t2),(t3,t4)),(t5,t6));")


# ===== T-393: an ambiguity token is a state set, not a state =================
# Deriving a secondary's levels from the observed token STRINGS admitted "{01}"
# as a level of its own, one Hamming step from both "0" and "1" rather than
# matching either -- so a polymorphic cell cost a step no resolution of it
# needs, and every polymorphic cell multiplied the combination count.

test_that("Polymorphic secondary costs no more than its best resolution", {
  mat <- matrix(c(
    "1", "0",
    "1", "1",
    "1", "{01}",
    "1", "0",
    "0", "-",
    "1", "1"
  ), nrow = 6, byrow = TRUE, dimnames = list(paste0("t", 1:6), NULL))
  ds <- xss_dat(mat)
  h <- CharacterHierarchy("1" = 2L)
  tree <- xss_tree()

  # "{01}" is a subset of {"0", "1"}, so the block's length is the smallest any
  # concrete resolution attains -- never more.
  resolved <- vapply(c("0", "1"), function(state) {
    m <- mat
    m[3, 2] <- state
    TreeLength(tree, xss_dat(m), hierarchy = h, inapplicable = "xform")
  }, numeric(1))
  expect_equal(
    TreeLength(tree, ds, hierarchy = h, inapplicable = "xform"),
    min(resolved)
  )

  # The binary secondary has two levels, not three.
  expect_equal(RecodeHierarchy(ds, h)$sankoff_chars[[1]]$n_states, 3)
})


test_that("Polymorphic cells do not inflate the state space", {
  # Four binary secondaries = 2^4 + 1 = 17 states.  Reading a state space off
  # the token strings made each polymorphic cell a third level, giving 82 and a
  # spurious "> 32 states" warning (with the quadratic-per-node cost to match).
  mat <- matrix(c(
    "1", "0",    "0",    "0",    "0",
    "1", "1",    "1",    "1",    "1",
    "1", "{01}", "0",    "1",    "0",
    "1", "0",    "{01}", "0",    "1",
    "1", "1",    "0",    "{01}", "1",
    "1", "0",    "1",    "0",    "{01}",
    "0", "-",    "-",    "-",    "-"
  ), nrow = 7, byrow = TRUE, dimnames = list(paste0("t", 1:7), NULL))
  ds <- xss_dat(mat)
  h <- CharacterHierarchy("1" = 2:5)

  expect_silent(recoded <- RecodeHierarchy(ds, h))
  expect_equal(recoded$sankoff_chars[[1]]$n_states, 17)
})


test_that("Polymorphism narrows a multistate secondary without freeing it", {
  # With three levels available, "{01}" is neither resolved nor unconstrained:
  # t2 must not be allowed to take state "2" to match its sister t1.
  mat <- matrix(c(
    "1", "2",
    "1", "{01}",
    "1", "2",
    "1", "0",
    "0", "-",
    "1", "1"
  ), nrow = 6, byrow = TRUE, dimnames = list(paste0("t", 1:6), NULL))
  ds <- xss_dat(mat)
  h <- CharacterHierarchy("1" = 2L)
  tree <- xss_tree()

  resolved <- vapply(c("0", "1", "2"), function(state) {
    m <- mat
    m[2, 2] <- state
    TreeLength(tree, xss_dat(m), hierarchy = h, inapplicable = "xform")
  }, numeric(1))
  # Precondition: resolving to "2" is strictly cheaper here, so treating the
  # token as wholly unknown would be a measurable under-count.
  expect_lt(resolved[["2"]], min(resolved[c("0", "1")]))

  expect_equal(
    TreeLength(tree, ds, hierarchy = h, inapplicable = "xform"),
    min(resolved[c("0", "1")])
  )
})


# ===== T-394: a secondary with no observed state ============================
# `ValidateHierarchy()` runs before both taxon-subsetting sites, so dropping a
# taxon can leave a secondary all-gap/all-missing in a dataset that validated.
# A zero-width state space then admitted no state at any tip: every tip cost
# was infinite and the block's length `Inf`.

xss_degenerate <- function() {
  # Character 3 is resolved at s5 alone; scoring any tree over {s1..s4} drops
  # s5, leaving char 3 with nothing observed.
  matrix(c(
    "1", "0", "?",
    "1", "1", "?",
    "1", "0", "?",
    "0", "-", "-",
    "1", "1", "0"
  ), nrow = 5, byrow = TRUE, dimnames = list(paste0("s", 1:5), NULL))
}

test_that("Unobserved secondary leaves a finite length after taxon dropping", {
  ds <- xss_dat(xss_degenerate())
  h <- CharacterHierarchy("1" = 2:3)
  tree <- ape::read.tree(text = "((s1,s2),(s3,s4));")

  # Precondition: the dataset as supplied validates and recodes normally.
  expect_equal(RecodeHierarchy(ds, h)$sankoff_chars[[1]]$n_states, 3)

  expect_true(is.finite(
    TreeLength(tree, ds, hierarchy = h, inapplicable = "xform")))
})


test_that("Search over a subset with an unobserved secondary completes", {
  ds <- xss_dat(xss_degenerate())
  h <- CharacterHierarchy("1" = 2:3)
  tree <- ape::read.tree(text = "((s1,s2),(s3,s4));")

  res <- suppressWarnings(
    MaximizeParsimony(ds, tree = tree, hierarchy = h, inapplicable = "xform",
                      maxReplicates = 2L, verbosity = 0L))
  expect_true(is.finite(attr(res, "score")))
})


test_that("A pool with no finite length is reported, not branched on", {
  # `diff(range(c(Inf, Inf)))` is `NaN`, which aborted the reporting `if` with
  # "missing value where TRUE/FALSE needed" whatever produced the infinities.
  expect_warning(reported <- TreeSearch:::.ReportXformScore(c(Inf, Inf), 7),
                 "no finite x-transformation length")
  expect_equal(reported, 7)

  # A pool that is only partly infinite still reports the shortest finite
  # length (and, these two differing, warns about that too).
  expect_warning(
    expect_warning(mixed <- TreeSearch:::.ReportXformScore(c(Inf, 3, 5), 7),
                   "no finite x-transformation length"),
    "do not share a length")
  expect_equal(mixed, 3)
})


# ===== T-401: the report block must not score a contracted tree =============
# `MaximizeParsimony()` rescores its own output at the canonical rooting.  With
# `collapse = TRUE` that output may be polytomous, and `TreeLength()` scores the
# topology it is given: the length of a polytomy is that of its best resolution,
# which the Sankoff kernel does not compute.  Only the HSJ/XFORM no-op in
# `compute_collapsed_flags_aggressive()` keeps the returned trees binary today.

test_that("Contracted trees are not scored as though binary", {
  # The T-330 fixture, whose equal-weights arm genuinely collapses (10 -> 8
  # edges) while the hierarchy character supports the contracted clade.
  mat <- matrix(c(
    "0", "0",
    "0", "0",
    "1", "0",
    "1", "1",
    "1", "1",
    "1", "-"
  ), nrow = 6, byrow = TRUE, dimnames = list(paste0("t", 1:6), NULL))
  ds <- phangorn::phyDat(mat, type = "USER", levels = c("-", "0", "1"),
                         ambiguity = "?")
  h <- CharacterHierarchy("2" = integer(0))
  binary <- Preorder(RenumberTips(
    ape::read.tree(text = "(((t1,t2),(t3,(t4,t5))),t6);"), names(ds)))

  at <- attributes(ds)
  collapsed <- TreeSearch:::ts_collapse_pool(
    list(binary[["edge"]]), at$contrast,
    matrix(unlist(ds, use.names = FALSE), nrow = length(ds), byrow = TRUE),
    as.integer(TreeSearch:::.NonHierarchyWeights(ds, h)), at$levels,
    list(min_steps = integer(0), concavity = Inf, xpiwe = FALSE,
         xpiwe_r = 0.5, xpiwe_max_f = 5.0, obs_count = integer(0),
         infoAmounts = NULL),
    NULL, NULL, NULL)
  polytomous <- Renumber(structure(
    list(edge = collapsed$trees[[1]],
         Nnode = max(collapsed$trees[[1]]) - length(ds),
         tip.label = names(ds)),
    class = "phylo"))
  expect_lt(dim(polytomous[["edge"]])[[1]], dim(binary[["edge"]])[[1]])

  reference <- TreeLength(binary, ds, hierarchy = h, inapplicable = "xform")

  # Precondition: handing the contracted tree straight to `TreeLength()` does
  # not reproduce that length.  Tolerant of an error, because the kernel-level
  # boundary check on non-binary edge matrices (T-400) turns the silent form of
  # this into a loud one.
  direct <- tryCatch(
    TreeLength(structure(list(polytomous), class = "multiPhylo"), ds,
               hierarchy = h, inapplicable = "xform"),
    error = function(e) NA_real_)
  expect_false(isTRUE(all.equal(direct, reference)))

  # The reporting path scores the binary pool the tree was contracted from.
  expect_equal(
    TreeSearch:::.XformPoolScore(list(polytomous), list(binary), ds, h, -1),
    reference)
})
