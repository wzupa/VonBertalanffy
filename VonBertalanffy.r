############################################################################
#   Von Bertalanffy: R code to estinate the Von Bertalanffy growth curve   #
#   Authors: W. Zupa                                                       #
#   Fondazione Coispa                                                      #
#   If you have any comments or suggestions please contact the following   #
#   e-mail address: zupa@fondazionecoispa.org                              #
#   June 2023                                                              #
############################################################################

# combined sexes ('C') are composed aggregating all the available sexes.
# The analysis by sex is conducted aggregating undetermined sex ('I') with
# each other (e.g. M = M + I; F = F + I).

library(ggplot2)

# set the working directory where to save results
setwd("D:\\OneDrive - Coispa Tecnologia & Ricerca S.C.A.R.L\\BLACK SEA\\Von Bertalanffy")

# load the data fron csv file (use ';' as value separator)
data <- read.table("data_example.csv", sep = ";", header = TRUE)

#---  DON'T MODIFY THE FOLLOWING CODE ---#

wd <- getwd()
data$SEX <- toupper(data$SEX)
data_C <- data[!is.na(data$AGE) &
  as.character(data$AGE) != " " &
  as.character(data$AGE) != "" &
  as.character(data$AGE) != "UR" &
  as.character(data$AGE) != "NR" &
  as.character(data$AGE) != "-1" &
  data$SEX %in% c("F", "M", "N"), ]
data_C$SEX <- "C"

data <- data[!is.na(data$AGE) &
  as.character(data$AGE) != " " &
  as.character(data$AGE) != "" &
  as.character(data$AGE) != "UR" &
  as.character(data$AGE) != "NR" &
  as.character(data$AGE) != "-1" &
  data$SEX %in% c("F", "M", "I", "N"), ]

data <- rbind(data, data_C)

data$AGE <- as.numeric(data$AGE)
data <- data[!is.na(data$AGE), ]

sp <- unique(data$SPECIES)
if (length(sp) > 1) {
  stop("Data for more than one species are not allowed")
}

vbf3 <- function(t, Linf, k, t0) {
  Linf * (1 - exp(-k * (t - t0)))
}

if (nrow(data) > 1) {
  sexes <- unique(data$SEX)
  s <- 2
  params <- list()
  for (s in 1:length(sexes)) {
    df <- data[data$SEX %in% c(sexes[s], "I"), ]
    if (nrow(df) > 0) {
      ss <- list()
      ss$Linf <- max(df$LENGTH) / 0.95
      ss$K <- mean(df[df$AGE == 1, "LENGTH"]) / ss$Linf
      ss$t0 <- -0.5

      error <- try(mod3 <- nls(LENGTH ~ vbf3(AGE, Linf, k, t0), start = c(Linf = abs(ss$Linf), k = abs(ss$K), t0 = ss$t0), data = df, algorithm = "default"), silent = TRUE)
      if (!is(error, "try-error")) {
        mod3 <- nls(LENGTH ~ vbf3(AGE, Linf, k, t0), start = c(Linf = abs(ss$Linf), k = abs(ss$K), t0 = ss$t0), data = df, algorithm = "default")
        sum <- summary(mod3)
        sum
        AIC(mod3)
        Linf <- sum$parameters[1]
        k <- sum$parameters[2]
        t0 <- sum$parameters[3]
        params[[s]] <- data.frame(species = sp, sex = sexes[s], Linf = as.numeric(round(Linf, 4)), k = as.numeric(round(k, 4)), t0 = as.numeric(round(t0, 4)), notes = "")
        pred <- data.frame(AGE = seq(0, max(df$AGE), 0.1))
        pred$LENGTH <- predict(mod3, newdata = data.frame(AGE = pred$AGE))
        p <- ggplot() +
          geom_point(data = df, aes(x = AGE, y = LENGTH), colour = "red", stat = "identity") +
          geom_line(data = pred, aes(x = AGE, y = LENGTH), colour = "blue", size = 1) +
          ggtitle(paste("VBGC", sexes[s], sp, sep = " - ")) +
          labs(col = "Year") +
          xlab("Age") +
          theme(legend.position = "none") +
          ylab(paste("Length (mm)"))
        print(p)
        ggsave(file.path(wd, "output", "VBGC", paste0("VBGC_", sp, "_", sexes[s], ".jpg")), dpi = 300, width = 9, height = 6, units = "in")
      } else {
        params[[s]] <- data.frame(species = sp, sex = sexes[s], Linf = NA, k = NA, t0 = NA, notes = "model not converged")
        p <- ggplot() +
          geom_point(data = df, aes(x = AGE, y = LENGTH), colour = "red", stat = "identity") +
          ggtitle(paste("VBGC", sexes[s], sp, sep = " - ")) +
          labs(col = "Year") +
          xlab("Age") +
          theme(legend.position = "none") +
          ylab(paste("Length (mm)"))
        print(p)
        ggsave(file.path(wd, "output", "VBGC", paste0("VBGC_", sp, "_", sexes[s], ".jpg")), dpi = 300, width = 9, height = 6, units = "in")
      }
    } else {
      cat(paste0("Not age data for '", sexes[s], "' sex"))
    }
  }
  df.params <- do.call(rbind, params)
  write.table(df.params, file.path(wd, "output", "VBGC", paste0("VBGC_", sp, "_summary_table.csv")), sep = ";", row.names = FALSE)
  print(df.params)
} else {
  message("Not enougth data to plot VBGC")
}
