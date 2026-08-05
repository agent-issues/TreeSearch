#' Recode hierarchical characters as step-matrix characters
#'
#' Implements the x-transformation recoding of
#' \insertCite{Goloboff2021;textual}{TreeSearch}.
#' Each hierarchy block (one controlling primary character plus \eqn{n}
#' secondary characters) is combined into a single step-matrix character
#' with \eqn{\prod \max(k_i, 1) + 1} states and an asymmetric cost matrix.
#'
#' @details
#' ## State encoding
#'
#' State 0 represents "primary absent".
#' States \eqn{1 \ldots \prod \max(k_i, 1)} represent all possible combinations
#' of secondary character states (where \eqn{k_i} is the number of informative
#' states of secondary character \eqn{i}; a secondary with none contributes a
#' single unobserved state, as below).
#'
#' The informative levels of a secondary character are read from the dataset's
#' `contrast` matrix, not from the token strings it carries.  An ambiguity token
#' such as `"{01}"` therefore denotes the *set* of states it contrasts against,
#' as it does under `inapplicable = "hsj"`, rather than becoming a level of its
#' own; a token that admits every applicable state (`"?"`, or an ambiguity
#' spanning them all) shows that no state is established and so contributes
#' none.  A secondary whose levels are all unobserved -- reachable by dropping
#' taxa from a dataset that validated -- is carried as a single unobserved
#' level: it adds nothing to any tree's length, but still counts towards the
#' block's gain cost, which is a property of the hierarchy rather than of the
#' taxa sampled.  A single secondary may take at most 31 levels.
#'
#' ## Cost matrix
#'
#' - **Absent → present (gain):** cost = \eqn{n + 1}, where \eqn{n} is the
#'   number of secondary characters.
#' - **Present → absent (loss):** cost = 1.
#' - **Present → present:** Hamming distance (number of secondaries with
#'   different states).
#'
#' ## Rooting
#'
#' Because gain and loss cost differently, this cost matrix is **asymmetric**
#' whenever a block has at least one secondary character -- and the length of a
#' tree under an asymmetric step matrix depends on where the tree is rooted,
#' unlike ordinary parsimony.  The asymmetry is the point of the recoding (the
#' first gain of the controlling character pays for the secondaries it brings
#' into existence), so this is intrinsic rather than a defect.
#'
#' `TreeSearch` treats topologies as unrooted, so [`TreeLength()`] and
#' [`MaximizeParsimony()`] both evaluate x-transformation lengths at a canonical
#' rooting -- on the first taxon of `dataset` -- giving one length per topology
#' and making a reported score reproducible.  That length is an upper bound on
#' the rooting-free minimum, exceeding it by at most the total number of
#' secondary characters across blocks.  If a rooting is biologically meaningful
#' to you, score the tree yourself with the block's `cost_matrix` rather than
#' relying on the canonical value.
#'
#' @param dataset A [`phyDat`][phangorn::phyDat] object.
#' @param hierarchy A [`CharacterHierarchy`] object.
#'
#' @return A list with elements:
#' \describe{
#'   \item{`sankoff_chars`}{A list of per-block lists, each containing:
#'     \describe{
#'       \item{`n_states`}{Integer, number of states (absent + present combos).}
#'       \item{`cost_matrix`}{Numeric matrix (\code{n_states × n_states}),
#'         row-major: \code{cost_matrix[from, to]}.}
#'       \item{`tip_states`}{Integer vector (length \code{n_tip}, 0-based).
#'         0 = absent, 1..n_present = present combination,
#'         -1 = fully ambiguous (all states possible),
#'         -2 = present but unknown combination.}
#'       \item{`forced_root_state`}{Integer: -1 (unconstrained).}
#'       \item{`block_chars`}{Integer vector of original character indices
#'         (1-based) belonging to this block.}
#'       \item{`combo_grid`}{Integer matrix (\code{n_present × n_secondary}),
#'         row \code{i} giving the 1-based level index of each secondary for
#'         present-state \code{i + 1}.}
#'       \item{`tip_sec_known`}{Integer matrix (\code{n_tip × n_secondary}).
#'         For tips with \code{tip_states == -2}, column \code{s} holds a
#'         bit mask of the levels secondary \code{s} may take at that tip
#'         (bit \code{i - 1} set = level \code{i} admissible), or 0 where it
#'         is unconstrained; used to restrict the admissible states of a
#'         partially-known combination.}
#'     }
#'   }
#'   \item{`non_hierarchy_indices`}{Integer vector of original character
#'     indices (1-based) not in any hierarchy block.}
#' }
#'
#' @references
#' \insertAllCited{}
#' @family tree scoring
#' @seealso [CharacterHierarchy()], [MaximizeParsimony()]
#' @keywords internal
#' @export
RecodeHierarchy <- function(dataset, hierarchy) {
  ValidateHierarchy(hierarchy, dataset)

  idx <- attr(dataset, "index")
  allLevels <- attr(dataset, "allLevels")
  levels <- attr(dataset, "levels")
  contrast <- attr(dataset, "contrast")
  nChar <- length(idx)
  nTip <- length(dataset)

  # Original character matrix (taxon × char), as `contrast` row indices -- which
  # is what a secondary's state space must be read from.  Reading it off the
  # token strings instead made an ambiguity token such as "{01}" a level of its
  # own: one Hamming step from both "0" and "1" rather than matching either, and
  # one more factor in the combination count (T-393).  `tokenLevels` is the
  # R-side counterpart of `DataSet::token_states`, which the HSJ path reads.
  tokenMat <- do.call(rbind, lapply(dataset, function(x) x[idx]))
  applicable <- which(levels != "-")
  tokenLevels <- lapply(seq_len(nrow(contrast)), function(tk) {
    applicable[contrast[tk, applicable] > 0]
  })
  # A token admitting every applicable state establishes no state at all; this
  # is the role "?" played under the old string test, and an ambiguity spanning
  # the whole state space says exactly as much.
  tokenGeneric <- lengths(tokenLevels) == length(applicable)

  .RecodeBlock <- function(node) {
    ctrl <- node$controlling
    deps <- node$dependents

    if (length(node$children) > 0L) {
      stop("Nested hierarchies not yet supported in RecodeHierarchy(). ",
           "Block controlled by character ", ctrl, " has sub-hierarchies.")
    }

    # Informative levels for each secondary, as state indices into `levels`
    secLevels <- lapply(deps, function(d) {
      tokens <- unique(tokenMat[, d])
      sort(unique(unlist(tokenLevels[tokens[!tokenGeneric[tokens]]])))
    })
    secNStates <- vapply(secLevels, length, integer(1))
    if (any(secNStates > 31L)) {
      stop("Secondary character ", deps[which.max(secNStates)],
           " has more than 31 informative states; the x-transformation ",
           "cannot recode it.")
    }
    # A secondary with no informative level -- every tip gap or fully ambiguous
    # -- is carried as one unobserved level rather than dropped, so that a
    # present primary still has a state to take (a zero-width state space made
    # every tip cost infinite, T-394) and the block's gain cost still reflects
    # how many secondaries the primary controls, not how many the sampled taxa
    # happen to resolve.
    secNLevels <- pmax(secNStates, 1L)

    nPresent <- prod(secNLevels)
    nStates <- nPresent + 1L
    nSec <- length(deps)

    if (nStates > 32L) {
      warning(sprintf(
        paste0("Hierarchy block controlled by character %d produces %d states ",
               "(> 32). Large state spaces may be slow."),
        ctrl, nStates
      ))
    }

    # All present-state combinations (expand.grid: first dim varies fastest)
    if (nSec > 0L) {
      comboGrid <- as.matrix(expand.grid(
        lapply(secNLevels, seq_len)
      ))
    } else {
      # No secondaries: 2 states (absent + one present)
      comboGrid <- matrix(integer(0), nrow = 1L, ncol = 0L)
    }

    # --- Cost matrix ---
    gainCost <- nSec + 1L
    cm <- matrix(0, nStates, nStates)
    for (i in seq_len(nStates)) {
      for (j in seq_len(nStates)) {
        if (i == j) next
        if (i == 1L) {
          cm[i, j] <- gainCost  # absent → present
        } else if (j == 1L) {
          cm[i, j] <- 1         # present → absent
        } else {
          # Hamming distance between present combinations
          cm[i, j] <- sum(comboGrid[i - 1L, ] != comboGrid[j - 1L, ])
        }
      }
    }

    # --- Tip states ---
    # `tipSecKnown[t, s]` records, per tip and per secondary, a bit mask of the
    # levels that secondary may take at this tip (bit i - 1 = level i), or 0
    # where it is unconstrained. Only consulted when `tipStates[t] == -2`
    # (present, but not every secondary was resolvable): it lets the
    # admissible-state set be restricted to combinations consistent with
    # whatever the secondaries WERE observed to be, rather than freeing every
    # present state (T-379). A mask rather than a single level index because a
    # polymorphic token narrows a secondary without resolving it (T-393).
    tipStates <- integer(nTip)
    tipSecKnown <- matrix(0L, nrow = nTip, ncol = nSec)
    for (t in seq_len(nTip)) {
      pri <- allLevels[tokenMat[t, ctrl]]

      if (pri == "?") {
        tipStates[t] <- -1L  # fully ambiguous
        next
      }
      if (pri == "0" || pri == "-") {
        tipStates[t] <- 0L   # absent
        next
      }
      # Primary present: encode secondary combination
      if (nSec == 0L) {
        tipStates[t] <- 1L   # only present state
        next
      }

      secVals <- tokenMat[t, deps]
      anyUnknown <- FALSE
      levelIndices <- integer(nSec)
      secMasks <- integer(nSec)
      known <- logical(nSec)

      for (s in seq_len(nSec)) {
        if (tokenGeneric[[secVals[s]]]) {
          anyUnknown <- TRUE
          next
        }
        # Positions, within this secondary's levels, that its token admits.
        # `secLevels[[s]]` is the union over the non-generic tokens of this very
        # column, so a non-generic token's states are all levels of it and the
        # match cannot fail.
        pos <- match(tokenLevels[[secVals[s]]], secLevels[[s]])
        if (length(pos) == 1L) {
          levelIndices[s] <- pos
          known[s] <- TRUE
          next
        }
        anyUnknown <- TRUE
        # Admitting every level (or none of them) constrains nothing, and 0
        # says so more cheaply than the equivalent full mask.
        if (length(pos) > 0L && length(pos) < secNStates[[s]]) {
          secMasks[s] <- sum(bitwShiftL(1L, pos - 1L))
        }
      }

      if (anyUnknown) {
        tipStates[t] <- -2L  # present, one or more secondaries unresolved
        secMasks[known] <- bitwShiftL(1L, levelIndices[known] - 1L)
        tipSecKnown[t, ] <- secMasks
        next
      }

      # Mixed-radix encoding (first dim varies fastest, matching expand.grid)
      rowIdx <- 1L
      multiplier <- 1L
      for (s in seq_len(nSec)) {
        rowIdx <- rowIdx + (levelIndices[s] - 1L) * multiplier
        multiplier <- multiplier * secNLevels[s]
      }
      tipStates[t] <- rowIdx  # 1-based present state = Sankoff state index
    }

    list(
      n_states = nStates,
      cost_matrix = cm,
      tip_states = tipStates,
      forced_root_state = -1L,
      block_chars = c(ctrl, deps),
      # 1-based level index of each secondary for each present-state combo
      # (row i = state i + 1), used to resolve `tip_states == -2` sentinels.
      combo_grid = comboGrid,
      tip_sec_known = tipSecKnown
    )
  }

  blocks <- lapply(hierarchy, .RecodeBlock)

  hChars <- HierarchyChars(hierarchy)
  nonH <- setdiff(seq_len(nChar), hChars)

  list(
    sankoff_chars = blocks,
    non_hierarchy_indices = nonH
  )
}
