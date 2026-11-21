# Synthetic responses in the format produced by make_collector() (three-point
# elicitation on a 0-10 scale) for twelve fictitious members; used in examples.
set.seed(2026)
members <- c("Alba", "Bruno", "Carla", "Dario", "Elena", "Fabio", "Gaia", "Luca", "Marta", "Nico", "Olga", "Paolo")
grp <- c(1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3)
rows <- list()
for (i in seq_along(members)[-c(4, 11)]) for (j in seq_along(members)[-i]) {
  if (runif(1) > ifelse(grp[i] == grp[j], 0.85, 0.2)) next
  best <- round(min(10, max(0, rnorm(1, ifelse(grp[i] == grp[j], 7, 3), 1.5))) * 2) / 2
  u <- round(runif(2, 0, 2.5) * 2) / 2
  rows[[length(rows) + 1]] <- data.frame(group = "demo", timestamp = "2026-11-12T10:00:00Z",
    rater = members[i], rated = members[j], low = max(0, best - u[1]), best = best, high = min(10, best + u[2]))
}
d <- do.call(rbind, rows)
write.csv(d, "inst/extdata/three_point_example.csv", row.names = FALSE)
