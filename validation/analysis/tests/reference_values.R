# Reference values for test_stats.py, from metafor and pROC.
#   Rscript reference_values.R > reference_values.json
suppressMessages({library(metafor); library(pROC)})
y <- c(0.90, 0.95, 0.86, 0.80, 0.93, 0.71)
v <- c(0.0012, 0.0008, 0.0020, 0.0030, 0.0010, 0.0045)
# Tight convergence: the default threshold (1e-5) stops early
reml <- rma(y, v, method = "REML", test = "knha",
            control = list(threshold = 1e-12))
dl <- rma(y, v, method = "DL")
# Paired scores with ties
cases_a <- c(2.1, 3.4, 0.5, 1.9, 2.2, 4.0, 1.1, 2.2, 3.3, 0.9, 2.8, 1.5)
ctrl_a  <- c(-1.0, 0.2, 0.5, -0.3, 1.2, -2.0, 0.0, 0.8, 1.9, -0.5)
cases_b <- c(0.10, 0.30, 0.02, 0.12, 0.05, 0.40, 0.08, 0.20, 0.15, 0.01, 0.30, 0.09)
ctrl_b  <- c(0.00, 0.02, 0.05, 0.01, 0.10, 0.00, 0.00, 0.03, 0.20, 0.02)
resp <- c(rep(1, 12), rep(0, 10))
ra <- roc(resp, c(cases_a, ctrl_a), direction = "<", quiet = TRUE)
rb <- roc(resp, c(cases_b, ctrl_b), direction = "<", quiet = TRUE)
t <- roc.test(ra, rb, method = "delong", paired = TRUE)
num <- function(x) sprintf("%.12g", as.numeric(x))
cat("{\n")
cat(sprintf('"meta_y": [%s],\n', paste(num(y), collapse = ", ")))
cat(sprintf('"meta_v": [%s],\n', paste(num(v), collapse = ", ")))
cat(sprintf('"reml_tau2": %s, "reml_est": %s, "reml_se": %s, "reml_ci": [%s, %s],\n',
            num(reml$tau2), num(reml$b), num(reml$se), num(reml$ci.lb), num(reml$ci.ub)))
cat(sprintf('"dl_tau2": %s, "dl_est": %s, "dl_se": %s, "dl_ci": [%s, %s], "dl_i2": %s, "Q": %s,\n',
            num(dl$tau2), num(dl$b), num(dl$se), num(dl$ci.lb), num(dl$ci.ub),
            num(dl$I2 / 100), num(dl$QE)))
cat(sprintf('"cases_a": [%s], "ctrl_a": [%s], "cases_b": [%s], "ctrl_b": [%s],\n',
            paste(num(cases_a), collapse = ", "), paste(num(ctrl_a), collapse = ", "),
            paste(num(cases_b), collapse = ", "), paste(num(ctrl_b), collapse = ", ")))
cat(sprintf('"auc_a": %s, "auc_b": %s, "var_a": %s, "delong_z": %s, "delong_p": %s\n',
            num(auc(ra)), num(auc(rb)), num(var(ra, method = "delong")),
            num(t$statistic), num(t$p.value)))
cat("}\n")
