# form[[c(3, 1)]]' Table one summaries
#'
#' Summarizes baseline trial results by treatment
#' @param data dataframe
#' @param form formula y ~ x1 + x2
#' @param ... extra parameters passed through to speciality functions
#' @return a dataframe
#' @examples
#' table1(dat2, form = arm ~ sex + age, annot = FALSE)

#' @export
table1 <- function(form, data, ...) {
  UseMethod("table1")
}

row_name <- function(x, nm, ...) {
  UseMethod("row_name")
}

row_name.character <- function(x, nm) {
  categs <- unique(na.omit(as.character(x)))
  nms <- cbind(vars = c(nm, categs), code = c(1, rep(2, length(categs))))
  # browser()
  return(as.data.frame(nms))
}
row_name.factor <- row_name.character

row_name.logical <- row_name.character

row_name.numeric <- function(x, nm, ...) {
  return(as.data.frame(cbind(vars = nm, code = 3)))
}

row_summary <- function(x, yy) {
  UseMethod("row_summary")
}

row_summary.character <- function(x, yy) {
  df <- data.frame(x = x, y = yy) |> na.omit()
  t1 <- df |>
    tabyl(x, y, show_missing_levels = FALSE) |>
    adorn_percentages("col") |>
    adorn_pct_formatting(digits = 0) |>
    adorn_ns(position = "front") |>
    select(-x)
  browser()
  return(rbind("", t1))
}

row_summary.factor <- row_summary.character
row_summary.logical <- row_summary.character

row_summary.numeric <- function(x, yy) {
  sp <- split(x, yy)
  mm <- sp |>
    map_vec(mean, na.rm = TRUE) |>
    round(2)
  ss <- sp |>
    map_vec(sd, na.rm = TRUE) |>
    round(2) |>
    paste0("(", x = _, ")")
  out <- paste(mm, ss)
  names(out) <- names(sp)
  return(out)
}

row_pv <- function(x, yy) {
  UseMethod("row_pv")
}

row_pv.character <- function(x, yy) {
  tab <- data.frame(x = x, y = yy) |>
    na.omit() |>
    tabyl(x, y, show_missing_levels = FALSE)
  if (!(nrow(tab) >= 2 & ncol(tab) >= 3)) {
    pv <- NA
  } else {
    pv <- janitor::fisher.test(tab, simulate.p.value = TRUE)$p.value |>
      round(4)
  }
  return(c(pv, rep("", nrow(tab))))
}

row_pv.factor <- row_pv.character
row_pv.logical <- row_pv.character
row_pv.factor <- row_pv.character
row_pv.logical <- row_pv.character

row_pv.numeric <- function(x, yy) {
  categs <- unique(na.omit(yy))
  if (!(length(categs) > 1)) {
    return(NA)
  }
  df <- data.frame(x = x, y = yy)
  pv <- tidy(anova(lm(x ~ y, data = df)))$p.value[1] |>
    round(4)
  return(pv)
}

build <- function(indep, dep) {
  left <- indep |>
    imap(row_name) |>
    bind_rows()
  right <- indep |>
    map(row_pv, yy = dep[[1]]) |>
    unlist() |>
    enframe(name = NULL) |>
    setNames("p.value")
  mid <- indep |>
    map(row_summary, yy = dep[[1]]) |>
    bind_rows()
  names(mid) <- paste0(names(mid), " (N = ", table(dep), ")")

  tab <- bind_cols(left, mid, right)
  return(tab)
}
#' @export
#' @describeIn table1 interprets formula and yields publication tables
table1.formula <- function(form, data, ...) {
  args <- list(...)
  for (i in 1:length(args)) {
    assign(x = names(args)[i], value = args[[i]])
  }
  vars <- all.vars(form)
  y_var <- deparse(form[[2]])
  g_bar <- deparse(form[[c(3, 1)]])
  g_var <- NULL
  if (g_bar == "|") {
    x_vars <- all.vars(form[[c(3, 2)]])
    g_var <- all.vars(form[[c(3, 3)]])
    group <- data[g_var]
  } else {
    x_vars <- all.vars(form)[-1]
  }

  if (!is.null(g_var)) {
    dd <- split(data[x_vars], data[g_var])
    yy <- split(data[y_var], data[g_var])
    tab3 <- map2(dd, yy, ~ build(indep = .x, dep = .y))
    tab4 <- bind_rows(tab3)
    new <- data.frame(matrix(NA, nrow = length(tab3), ncol = ncol(tab4)))
    names(new) <- names(tab4)
    new$vars <- names(tab3)
    new$code <- 4
    rr <- cumsum(c(1, rep(nrow(tab3[[1]]) + 1, length(tab3) - 1)))
    tab5 <- insertRows(tab4, rr, new, rcurr = F)
  } else {
    tab5 <- build(indep = data[x_vars], dep = data[y_var])
  }
  stripes <- map(sort(unique(tab5$code)), function(x) {
    which(tab5$code == x)
  })
  myfcn <- function(x, i, theme = theme) {
    x <- x |> row_spec(stripes[[i]],
      color = theme$foreground[i],
      background = theme$background[i]
    )
  }
  if (layout == "console") {
    return(tab5[-2])
  } else if (layout == "latex") {
    names(tab5) <- gsub("\\(N", "\\\\(N", names(tab5))
    tab6 <- tab5 |> dplyr::mutate(
      vars = ifelse(code == 2, gsub("^", "\\\\quad ", vars), vars)
    )

    kk <- kbl(tab6[-2], "latex",
      booktabs = T, linesep = "",
      escape = F, digits = digits
    )
    tablatex <- reduce(1:length(stripes), ~ myfcn(.x, .y, theme = theme), .init = kk)
    write(tablatex, paste0("./tables/", fname, ".tex"))
    system(paste0("sh ~/shr/figurize.sh ./tables/", fname, ".tex"))
  } else if (layout == "html") {
    # browser()
    kk <- kbl(tab5[-2], "html",
      escape = F, digits = digits
    )
    tabhtml <- reduce(1:length(stripes), ~ myfcn(.x, .y, theme = theme), .init = kk)
  }
}

