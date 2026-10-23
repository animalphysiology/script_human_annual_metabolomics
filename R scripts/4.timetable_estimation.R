library(ggplot2)
library(patchwork)
library(RColorBrewer)
library(tidyverse)


# load data
load("data/Human_annual_metabolomics.RData")

source("code/supporting.functions.R")


seasonal.metabolites <- readRDS("results/human_amide_lmm_cosinor_annual.rds")
rownames(seasonal.metabolites) <- seasonal.metabolites$unique_id


################################################################################

LOO_lmm_cosinor <- readRDS("results/LOO_human_amide_lmm_cosinor_annual.rds")
LOO_lmm_cosinor <- lapply(LOO_lmm_cosinor, function(df) {
  df$unique_id <- row.names(df)
  df <- merge(human_amide_BGI[, 1:20], df, by = "unique_id")
  df
})

LOO_lmm_cosinor |>
  lapply(function(df) {
    subset(df, q.value < 0.05) |> nrow()
  }) |>
  unlist() |>
  min()

# 243

topN <- c(seq(10, 200, 10))

error <- list()

for (n in topN) {
  print(n)
  timetable.pred <- do.call(rbind, lapply(LETTERS[1:11], function(subject.out) {
    SOMs_sub <- subset(LOO_lmm_cosinor[[subject.out]], q.value < .05)
    SOMs_sub <- SOMs_sub[order(-log10(SOMs_sub$q.value), decreasing = TRUE), ]
    SOMs <- head(SOMs_sub, n)

    # divide data to test (1 subject) and train (10 subjects)
    test.metadata <- subset(sample_metadata, Subject == subject.out) |> as.data.frame()
    train.metadata <- subset(
      sample_metadata, Subject != subject.out
    ) |> as.data.frame()

    train <- log2(human_amide_BGI[, train.metadata$Sample])
    test <- log2(human_amide_BGI[, test.metadata$Sample])

    rownames(train) <- rownames(test) <- human_amide_BGI$unique_id

    test <- sweep(
      sweep(test, 1, rowMeans(train, na.rm = TRUE)),
      1, apply(train, 1, sd, na.rm = TRUE), "/"
    )

    for (month in test.metadata$Sample) {
      y <- test[SOMs$unique_id, month]
      fit <- glm(y ~ sin(phase * w) + cos(phase * w), data = SOMs)
      b <- coef(fit)

      p <- c(
        b,
        RS = 1 - deviance(fit) / deviance(glm(y ~ 1)),
        amplitude = amplitude(b[3], b[2]) |> as.numeric(),
        predicted = phase(b[3], b[2]) |> as.numeric(),
        actual = subset(test.metadata, Sample == month)$day_of_year %% period
      )

      p["error"] <- circular.diff(p["predicted"], p["actual"])
      test.metadata[test.metadata$Sample == month, names(p)] <- p
    }

    test.metadata
  }))

  error[[paste0("top", n)]] <- timetable.pred
}

################################################################################
colors <- colorRampPalette(
  colors = c("#0050FF", "grey50", "#FF5000")
)(length(topN))


timetable.errors <- do.call(rbind, error)
timetable.errors$top <-
  gsub("top([0-9]*)\\..*", "\\1", rownames(timetable.errors)) |>
  as.numeric() |>
  factor()


plt_RS_timetable_top <-
  group_by(timetable.errors, top) %>%
  summarise(
    "RS" = rsq(actual, actual + error)
  ) %>%
  ggplot() +
  geom_vline(xintercept = c(30), linetype = "dotted") +
  geom_point(aes(x = top |> as.character() |> as.numeric(), y = RS, color = top), size = 3) +
  scale_color_manual(values = colors) +
  scale_x_continuous(
    name = "Top N metabolites", breaks = c(seq(20, max(topN), 20)),
    limits = c(0, max(topN) + 10), expand = c(0, 0)
  ) +
  scale_y_continuous(
    name = expression("R"^2 * " of actual vs predicted sampling date"),
    limits = c(0.75, 1), expand = c(0, 0)
  ) +
  theme_annual(x.angle = 45) +
  theme(legend.position = "none")

plt_RS_timetable_top

# fig 3A
pdf("figures/timetable/timetable_q-value_top-SOMs_RS.pdf", height = 3, width = 4.5)
plt_RS_timetable_top
dev.off()


################################################################################
# top 30 was the best

timetable.pred <- error$top30

# coefficient of determination
with(
  timetable.pred,
  rsq(actual, actual + error)
)

# for significance testing, Dec, Jan, Feb data was jointed as month 0, 13, 14
# this make chronologically continuous date value

bs <- rbind(
  subset(timetable.pred, Month == 12),
  subset(timetable.pred, Month == 1),
  subset(timetable.pred, Month == 2)
)
bs$actual <- ifelse(bs$Month == 12, bs$actual - 365, bs$actual + 365)
bs$Month <- ifelse(bs$Month == 12, bs$Month - 12, bs$Month + 12)

bs <- rbind(timetable.pred, bs)
bs <- bs[order(bs$Subject), ]

lapply(0:11, function(n) {
  t.test(
    with(subset(bs, Month == n), actual + error),
    with(subset(bs, Month == (n + 1)), actual + error),
    paired = T
  )[["p.value"]]
}) |>
  p.adjust("holm") |>
  signif(2)

# Dec vs Jan, Jan vs Feb, ... May vs Jun
# 0.12000 0.12000 0.65000 0.00025 0.09100 0.12000
# Jun vs Jul, Jul vs Aug, ... Nov vs Dec
# 0.28000 0.00220 0.28000 0.28000 0.65000 0.65000

lapply(0:11, function(n) {
  t.test(
    with(subset(bs, Month == n), actual + error),
    with(subset(bs, Month == (n + 2)), actual + error),
    paired = T
  )[["p.value"]]
}) |>
  p.adjust("holm") |>
  signif(2)

# Dec vs Feb, Jan vs Mar, ... May vs Jul
# 6.1e-03 1.4e-03 1.8e-04 5.3e-06 1.1e-03 9.6e-04
# Jun vs Aug, Aug vs Oct, ... Nov vs Jan
# 2.7e-04 9.6e-04 1.1e-03 5.5e-02 3.1e-02 6.3e-03

lapply(0:11, function(n) {
  t.test(
    with(subset(bs, Month == n), actual + error),
    with(subset(bs, Month == (n + 3)), actual + error),
    paired = T
  )[["p.value"]]
}) |>
  p.adjust("holm") |>
  signif(2)

# all pairs were siginificant <0.05
# [1] 1.6e-04 3.8e-05 7.4e-05 1.3e-06 2.4e-06 1.3e-06 3.8e-05 3.8e-05 1.3e-02 4.4e-03 7.5e-05
# [12] 1.3e-03


pdf("figures/timetable/actual-estimated_fixed_top30_qval_ranked.pdf", width = 5, height = 5)

ggplot(data = subset(timetable.pred)) +
  geom_point(aes(
    x = actual,
    y = predicted,
    color = Subject,
    fill = Subject,
    shape = Subject
  ), size = 2) +
  scale_Subject() +
  scale_x_continuous(
    name = "actual sampling day", expand = c(0.01, 0.01),
    limits = c(0, period),
    sec.axis = sec_axis(
      ~ (. + 10) / period * 12,
      breaks = 1:12,
      labels = month.abb,
      name = "actual month"
    )
  ) +
  scale_y_continuous(
    name = "predicted sampling day", expand = c(0.01, 0.01),
    limits = c(0, period),
    sec.axis = sec_axis(
      ~ (. + 10) / period * 12,
      breaks = 1:12,
      labels = month.abb,
      name = "predicted month"
    )
  ) +
  geom_abline(intercept = 0, slope = 1, color = "grey20") +
  theme_annual() +
  coord_fixed()
dev.off()




################################################################################
# visualize timetable prediction
n <- 30
plots <- lapply(LETTERS[1:11], function(subject.out) {
  plots <- list()
  SOMs_sub <- subset(LOO_lmm_cosinor[[subject.out]], q.value < .05)
  SOMs_sub <- SOMs_sub[order(rank.function(SOMs_sub), decreasing = TRUE), ]
  SOMs <- head(SOMs_sub, n)

  test.metadata <- subset(sample_metadata, Subject == subject.out) |> as.data.frame()
  train.metadata <- subset(
    sample_metadata,
    Subject != subject.out # &Sex == ifelse(subject.out %in% LETTERS[1:6], "male", "female")
  ) |> as.data.frame()

  train <- log2(human_amide_BGI[, train.metadata$Sample])
  test <- log2(human_amide_BGI[, test.metadata$Sample])

  rownames(train) <- rownames(test) <- human_amide_BGI$unique_id

  test <- sweep(
    sweep(test, 1, rowMeans(train, na.rm = TRUE)),
    1, apply(train, 1, sd, na.rm = TRUE), "/"
  )

  for (month in test.metadata$Sample) {
    y <- test[SOMs$unique_id, month]
    fit <- glm(y ~ sin(phase * w) + cos(phase * w), data = SOMs)
    b <- coef(fit)

    p <- c(
      b,
      RS = 1 - deviance(fit) / deviance(glm(y ~ 1)),
      amplitude = amplitude(b[3], b[2]) |> as.numeric(),
      predicted = phase(b[3], b[2]) |> as.numeric(),
      actual = test.metadata$day_of_year[test.metadata$Sample == month] %% period
    )

    p["error"] <- circular.diff(p["predicted"], p["actual"])
    ylim <- c(min(y) |> floor(), max(y) |> ceiling())

    plots[[month]] <-
      ggplot(data = data.frame(y = y, x = SOMs$phase)) +
      geom_point(aes(x = x, y = y, color = x)) +
      geom_vline(xintercept = p["actual"], linetype = 2, color = "red") +
      geom_vline(xintercept = p["predicted"], linetype = 2, color = "black") +
      scale_color_gradientn(colors = circular.phase.color(1:365)) +
      stat_function(
        fun = eval(substitute(
          function(x) {
            A * cos((x - P) * W) + B
          },
          list(
            A = p["amplitude"],
            P = p["predicted"],
            W = w,
            B = p["(Intercept)"]
          )
        )),
        color = "black", linewidth = 0.8, n = period
      ) +
      scale_x_continuous(
        name = "reference peak phase of metabolites",
        breaks = month.label %% period,
        labels = (month.label %% period) |> names(),
        limits = c(0, period), expand = 0
      ) +
      scale_y_continuous(name = "relative abundance", limits = ylim, expand = 0) +
      theme_annual(x.angle = 45) +
      theme(legend.position = "n") +
      labs(title = paste(
        month.name[as.numeric(sub(".*_", "", month), "")],
        # "sample from subject", sub("_.*", "", month),
        ""
      ))
  }
  plots
})

pdf("figures/timetable/example_A_top30_qval_ranked.pdf", width = 12, height = 8)
wrap_plots(plots[[3]], ncol = 4, nrow = 3, byrow = F)
dev.off()
