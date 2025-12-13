source("../t2.R")

p_load(
  berryFunctions, atable, purrr, kableExtra, tibble,
  janitor, broom, palmerpenguins, dplyr
)
options(knitr.kable.NA = "")
p1 <- sample_n(penguins, 100) |>
  dplyr::select(
    species, flipper_length_mm, sex,
    body_mass_g, bill_length_mm, island
  ) |>
  dplyr::mutate(flp = flipper_length_mm > 197)
p1[100, "island"] <- NA

theme_npg <- list(
  foreground = c("black", "black", "black", "black"),
  background = c(
    "#f0efd4",
    "#e8e6bc",
    "#f0efd4",
    "#f0efd4"
  )
)
theme_nejm <- list(
  foreground = c("black", "black", "black", "black"),
  background = c("#fff7e9", "white", "#fff7e9", "white")
)

theme_green <- list(
  foreground = c("black", "black", "black", "black"),
  background = c("#99cfa8", "#d4f0dc", "#94ebad", "yellow")
)
theme_simple <- list(
  foreground = c("black", "black", "black", "black"),
  background = c("cyan", "blue", "green", "yellow")
)
theme_bw <- list(
  foreground = c("black", "black", "black", "black"),
  background = c("white", "white", "white", "white")
)
tab0 <- table1(sex ~ island,
  data = p1,
  theme = theme_green, layout = "console", fname = "ptab0", digits = 3, pvalue = FALSE
)
tab1 <- table1(sex ~ island + flp + body_mass_g + bill_length_mm,
  data = p1,
  theme = theme_green, layout = "latex", fname = "ptab1", digits = 3
)
tab1b <- table1(sex ~ island + flp + body_mass_g + bill_length_mm,
  data = p1,
  theme = theme_npg, layout = "latex", fname = "ptab1b", digits = 3
)
tab2 <- table1(sex ~ flp + body_mass_g + bill_length_mm | island,
  data = p1,
  theme = theme_npg, layout = "console", fname = "ptab2", digits = 3
)
tab3 <- table1(sex ~ flp + body_mass_g + bill_length_mm | island,
  data = p1,
  theme = theme_npg, layout = "latex", fname = "ptab3", digits = 3
)
tab4 <- table1(sex ~ flp + body_mass_g + bill_length_mm + island,
  data = p1,
  theme = theme_green, layout = "latex", fname = "ptab4", digits = 3, size=FALSE
)
tab4 <- table1(sex ~ flp + body_mass_g + bill_length_mm | island,
 data = p1,
 theme = theme_green, layout = "latex", fname = "ptab4b", digits = 3, size=FALSE
)
