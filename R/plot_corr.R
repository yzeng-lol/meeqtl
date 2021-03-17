################################################################################
############# correlation analysis related plotting functions ##################
##
##  corr_scatter: Scatterplot for x and y with correlation and lm results
##  corr_scatter_ci: Scatterplot for x and y with correlation and lm (with CI and Rugs)
##  corr_scatter_hex: Hexagonalfor x and y with correlation and lm
##  corr_volcano: Volcano plot for correlation coefficients and p values
##  corr_volcano_ylim: Volcano plot for correlation coefficients and p values with ylim
##  corr_dotmap


################################################################################
#' @title corr_scatter
#' scatterplot for x and y with correlation and lm information
#' recommend for light x and y
#'
#' @author yzeng
#' @param x,y A numeric vector
#' @param cor_m mehtod for correlation: "pearson", "kendall", "spearman"
#' @param x_lab The x axis lable
#' @param y_lab The y axis lable
#' @param title Name of output figure
#' @param fig_type The format of output figure, either "PNG/png" or "PDF/pdf"
#' @export
#'
#' @examples
#' x <- rnorm(30)
#' y <- rnorm(30)
#' corr_scatter(x, y, "pearson", "var_1", "var_2", "test_corr_scatter_plot", "PDF")

corr_scatter <- function(x, y, cor_m, x_lab, y_lab, title, fig_type){

    cor <- cor.test(x , y, method = cor_m)
    pval <- cor$p.value

    if(pval > 2.2e-16){
        pval <- formatC(pval, format = "e", digits = 2)
        sub_t <- paste(cor_m, ": r = ", round(cor$estimate, 2), "; ", "p_value = ", pval, sep = "")        # correlation coefficient
    } else {
        sub_t <- paste(cor_m, ": r = ", round(cor$estimate, 2), "; ", "p_value < 2.2e-16 ", sep = "")
    }

    ## png
    if (fig_type == "pdf" | fig_type == "PDF")
    {
        file_name <- paste(title, ".pdf", sep = "");
        pdf(file_name, width = 3, height = 3)
        par(mar = c(3, 3, 1, 1), mgp = c(2, 0.75, 0))
        plot(x, y, pch = 16, xlab = x_lab, ylab = y_lab, las = 1, main = sub_t, cex.main = 0.6, font.main = 3)
        abline(lm(y ~ x), lty = 2, col = "red")
        abline(a = 0, b = 1, lty = 2, col = "gray")
        dev.off()

    } else {

        file_name <- paste(title, ".png", sep = "");
        png(file_name, width = 3, height = 3, units = "in", res = 300)
        par(mar = c(3, 3, 1, 1), mgp = c(2, 0.75, 0))
        plot(x, y, pch = 16, xlab = x_lab, ylab = y_lab, las = 1, main = sub_t, cex.main = 0.6, font.main = 3)
        abline(lm(y ~ x), lty = 2, col = "red")
        abline(a = 0, b = 1, lty = 2, col = "gray")
        dev.off()
    }

}


################################################################################
#' @title corr_scatter_ci
#' scatterplot for x and y with correlation and lm (with CI and Rugs);
#' recommend for median numbers of  x and y
#'
#' @author yzeng
#' @param x,y A numeric vector
#' @param cor_m mehtod for correlation: "pearson", "kendall", "spearman"
#' @param x_lab The x axis lable
#' @param y_lab The y axis lable
#' @param file_name Name of output figure with format suffix
#' @import ggplot2
#' @improt ggpubr
#' @export
#'
#' @examples
#' x <- rnorm(100)
#' y <- rnorm(100)
#' corr_scatter_ci(x, y, "pearson", "var_1", "var_2", "test_corr_scatter_with_CI.pdf")

corr_scatter_ci <- function(x, y, cor_m,x_lab, y_lab, file_name)
{
    ## ggpubr::stat_cor can be an alternative way for adding r and p
    cor <- cor.test(x , y, method = cor_m)
    pval <- cor$p.value
    if(pval > 2.2e-16){
      pval <- formatC(pval, format = "e", digits = 2)
      sub_t <- paste(cor_m, ": r = ", round(cor$estimate, 2), "; ", "p_value = ", pval, sep = "")
    } else {
      sub_t <- paste(cor_m, ": r = ", round(cor$estimate, 2), "; ", "p_value < 2.2e-16 ", sep = "")
    }

    dat <- data.frame(x, y)
    g <- ggplot(dat, aes(x = x, y = y)) + geom_point(color='#2980B9', size = 2)
    g <- g + geom_smooth(method=lm, formula = y~x, color='#2C3E50') + theme_classic()
    g <- g + labs(x = x_lab, y = y_lab, title = sub_t) + theme(plot.title = element_text(size = 9))
    g <- g + geom_rug(color = "gray40", lwd = 0.35)
    ggsave(file_name, width = 3, height = 3)
}


################################################################################
#' @title corr_scatter_hex
#' Hexagonalfor x and y with correlation and lm
#' recommend for massive and compact x and y
#'
#' @author yzeng
#' @param x,y A numeric vector
#' @param cor_m mehtod for correlation: "pearson", "kendall", "spearman"
#' @param x_lab The x axis lable
#' @param y_lab The y axis lable
#' @param bin_cnt bins splict for x and y axis
#' @param file_name Name of output figure with format suffix
#' @import ggplot2
#' @improt ggpubr
#' @export
#'
#' @examples
#' x <- rnorm(10000)
#' y <- rnorm(10000)
#' corr_scatter_hex(x, y, "pearson", "var_1", "var_2", 50, "test_corr_scatter_hex.pdf")

corr_scatter_hex <- function(x, y, cor_m, x_lab, y_lab, bin_cnt, file_name)
{
  ## ggpubr::stat_cor can be an alternative way for adding r and p
  cor <- cor.test(x , y, method = cor_m)
  pval <- cor$p.value
  if(pval > 2.2e-16){
    pval <- formatC(pval, format = "e", digits = 2)
    sub_t <- paste(cor_m, ": r = ", round(cor$estimate, 2), "; ", "p_value = ", pval, sep = "")
  } else {
    sub_t <- paste(cor_m, ": r = ", round(cor$estimate, 2), "; ", "p_value < 2.2e-16 ", sep = "")
  }

  dat <- data.frame(x, y)
  g <- ggplot(dat, aes(x = x, y = y)) + geom_hex(bins = bin_cnt, show.legend = T)
  g <- g + geom_smooth(method=lm, formula = y~x, se = F, color='brown2', linetype="dashed", size = 0.75)
  g <- g + labs(x = x_lab, y = y_lab, title = sub_t) + theme_classic()
  ggsave(file_name, width = 3.6, height = 3)

}


################################################################################
#' @title corr_volcano
#' Volcano plot for correlation coefficients and p values
#'
#' @author yzeng
#' @param cor_r A numeric vector of correlation coefficients
#' @param cor_p A numeric vector of correlation p-values
#' @param name Name of output PNG figure
#' @param r_c The absuloate cutoff of correlation coefficient, default = 0.3
#' @param p_c The cor.test p-value cutoff, default = 0.05
#' @export
#'
#' @examples
#' cor_r <- runif(1000, -1, 1)
#' cor_p <- runif(1000, 0, 1)
#' corr_volcano(cor_r, cor_p, "test_corr_volcano_plot.png")

corr_volcano <- function(cor_r, cor_p, name, r_c = 0.3, p_c = 0.05)
{
  neg_r_c <- -r_c
  cor_p_log <- -log10(cor_p)
  ymax <- max(cor_p_log)
  xtxt_p = r_c + 0.2                 # x position for pos_correlated counts
  xtxt_n = -0.2 - r_c                # x position for neg_correlated counts
  ytxt_sig = 0.75*ymax               # y position for sig_cor count text
  ytxt_nosig = -log10(p_c)/2         # y position for nonsig_cor count test

  L <- length(cor_r)
  idx_u <- cor_r > r_c & cor_p < p_c
  idx_d  <- cor_r <  neg_r_c & cor_p < p_c
  idx_nu <- cor_r >= 0 & cor_r <= r_c & cor_p >= p_c
  idx_nd  <- cor_r <= 0 & cor_r >= neg_r_c & cor_p >= p_c

  ## Fisher-exact test
  t <- matrix(c(sum(idx_u, na.rm = T), sum(idx_d, na.rm = T),
                sum(idx_nu, na.rm = T), sum(idx_nd, na.rm = T)), 2, 2)
  tc <- fisher.test(t)
  pval <- as.numeric(tc$p.value)

  if(pval > 2.2e-16){
    pval <- formatC(pval, format = "e", digits = 2)
    sub_t <- paste("Fisher.test: p_value = ", pval, sep = "")
  } else {
    sub_t <- paste("Fisher.test: p_value < 2.2e-16")
  }

  ## color
  col_c <- rep("gray", L)
  col_c[idx_u] <- "red"
  col_c[idx_d] <- "blue"
  col_cnt <- table(col_c)

  png(file = name, width = 3, height = 3, units = "in", res = 600)
  par(mar = c(3, 3, 1, 1), mgp = c(2, 0.75, 0))

  plot(cor_r, cor_p_log, col = col_c, main = sub_t, cex.main = 0.6, las = 1,
       font.main = 3, xlab = "Correlaton Coefficient", ylab = "-Log10(p-value)")
  abline(h = -log10(p_c), lty = 2)
  segments(neg_r_c, -log10(p_c), x1 = neg_r_c, y1 = ymax, lty = 2)
  segments(r_c, -log10(p_c), x1 = r_c, y1 = ymax, lty = 2)
  segments(0, 0, x1 = 0, y1 = -log10(p_c), lty = 2)

  text(x = c(xtxt_n, xtxt_p), y = c(ytxt_sig, ytxt_sig), col = c("blue", "red"),
       labels = c(sum(idx_d, na.rm = T), sum(idx_u, na.rm = T)), cex = 0.75)

  text(x = c(xtxt_n, xtxt_p), y = c(ytxt_nosig, ytxt_nosig), col = c("gray40", "gray40"),
       labels = c(sum(idx_nd, na.rm = T), sum(idx_nu, na.rm = T)), cex = 0.75)

  dev.off()

}



################################################################################
#' @title corr_volcano_ylim
#' Volcano plot for correlation coefficients and p values. To deal with outliers
#' are with extreme p values, therefore if(y > ylim){y = ylim}
#'
#' @author yzeng
#' @param cor_r A numeric vector of correlation coefficients
#' @param cor_p A numeric vector of correlation p-values
#' @param y_lim The up limitation of the y
#' @param name Name of output PNG figure
#' @param r_c The absuloate cutoff of correlation coefficient, default = 0.3
#' @param p_c The cor.test p-value cutoff, default = 0.05
#' @export
#'
#' @examples
#' cor_r <- runif(1000, -1, 1)
#' cor_p <- runif(1000, 0, 1)
#' corr_volcano_ylim(cor_r, cor_p, 2, "test_corr_volcano_ylim_plot.png")

corr_volcano_ylim <- function(cor_r, cor_p, y_lim, name, r_c = 0.3, p_c = 0.05)
{
  neg_r_c <- -r_c
  cor_p_log <- -log10(cor_p)
  cor_p_log[cor_p_log > y_lim] <- y_lim       # if(y > ylim){y = ylim}
  ymax <- max(cor_p_log)

  xtxt_p = r_c + 0.2                 # x position for pos_correlated counts
  xtxt_n = -0.2 - r_c                # x position for neg_correlated counts
  ytxt_sig = 0.75*ymax               # y position for sig_cor count text
  ytxt_nosig = -log10(p_c)/2         # y position for nonsig_cor count test

  L <- length(cor_r)
  idx_u <- cor_r > r_c & cor_p < p_c
  idx_d  <- cor_r <  neg_r_c & cor_p < p_c
  idx_nu <- cor_r >= 0 & cor_r <= r_c & cor_p >= p_c
  idx_nd  <- cor_r <= 0 & cor_r >= neg_r_c & cor_p >= p_c

  ## Fisher-exact test
  t <- matrix(c(sum(idx_u, na.rm = T), sum(idx_d, na.rm = T),
                sum(idx_nu, na.rm = T), sum(idx_nd, na.rm = T)), 2, 2)
  tc <- fisher.test(t)
  pval <- as.numeric(tc$p.value)

  if(pval > 2.2e-16){
    pval <- formatC(pval, format = "e", digits = 2)
    sub_t <- paste("Fisher.test: p_value = ", pval, sep = "")
  } else {
    sub_t <- paste("Fisher.test: p_value < 2.2e-16")
  }

  ## color
  col_c <- rep("gray", L)
  col_c[idx_u] <- "red"
  col_c[idx_d] <- "blue"
  col_cnt <- table(col_c)

  png(file = name, width = 3, height = 3, units = "in", res = 600)
  par(mar = c(3, 3, 1, 1), mgp = c(2, 0.75, 0))
  plot(cor_r, cor_p_log, col = col_c, main = sub_t, cex.main = 0.6, las = 1,
       font.main = 3, xlab = "Correlaton Coefficient", ylab = "-Log10(p-value)")
  abline(h = -log10(p_c), lty = 2)
  segments(neg_r_c, -log10(p_c), x1 = neg_r_c, y1 = ymax, lty = 2)
  segments(r_c, -log10(p_c), x1 = r_c, y1 = ymax, lty = 2)
  segments(0, 0, x1 = 0, y1 = -log10(p_c), lty = 2)

  text(x = c(xtxt_n, xtxt_p), y = c(ytxt_sig, ytxt_sig), col = c("blue", "red"),
       labels = c(sum(idx_d, na.rm = T), sum(idx_u, na.rm = T)), cex = 0.75)
  text(x = c(xtxt_n, xtxt_p), y = c(ytxt_nosig, ytxt_nosig), col = c("gray40", "gray40"),
       labels = c(sum(idx_nd, na.rm = T), sum(idx_nu, na.rm = T)), cex = 0.75)
  dev.off()

}


################################################################################
#' @title corr_dotmap
#' Dotmap based on correlation coefficients and p values
#' @param cor_r A vector or matrix of correltaion coefficient
#' @param cor_p A vector or matirx of correlation test p values
#' @param r_cut,p_cut cutoffs for the significant correlations
#' @import corrplot
#' @export corr_dotmap
#' @examples
#'

corr_dotmap <- function(cor_r, cor_p, name){}
