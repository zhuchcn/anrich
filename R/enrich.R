#' @title enrichment test
#'
#' @description  Perform enrichment test on a concentration/expression matrix
#' of a hight through-put experiment study.
#'
#' A linear model was first fitted for each variable against the provided design
#' matrix. The statistics were then used to perform either the hypergeometric
#' test (Fisher's exact test) or Kolmogorov-Smirnov test to see whether a
#' particular category of variables are enriched.
#'
#' The function \code{\link{phyper}} is used for hypergeometric test for each
#' category in the parameter \strong{group}. For any given category \strong{g}
#' in \strong{group}, parameters for \code{\link{phyper}} are defined as
#' following. m is the total number of variables belong to category \strong{g}.
#' n is the total number of variables not belong to category \strong{g}. k is
#' the total number of variables that were affected in the study. x is the
#' number of variables belong to category \strong{g} and also affected.
#' Variables are defined as affected if the estimated coefficient is higher or
#' lower according to the parameter \strong{alternative}, with a p-value
#' smaller than the cutoff (p.cutoff).
#'
#' The function \code{\link{ks.test}} is used for Kolmogorov-Smirnov test for
#' each category in \strong{group}. The pvalues were first adjusted to one-tail
#' using \code{\link{pt}} according to the parameter \strong{alternative}. The
#' Kolmogorov-Smirnov test is then performed between the adjusted p-values of
#' variables that belong to the category \strong{g} and a sequence of simulated
#' values from 0 to 1.
#'
#' @param X The concentration/expression matrix. Rows are variables and columns
#' are observations.
#' @param design The design matrix often generated from
#' \code{\link{model.matrix}}. Number of rows must equal to the number of
#' columns of X.
#' @param group character. The categories of variables. Length must equal to the
#' number of rows of X.
#' @param coef character. The main effect of interest in the design matrix. Must
#' be one of the column names of the design matrix.
#' @param test character. Either hyper or ks.
#' @param p.cutoff numeric. The p-value cutoff to be used in hyper test. Default
#' is 1.
#'
#' @return A object of class either "enrichment.hyper" or "enrichment.ks" is
#' returned by this function. The returned object has the elements as described
#' below.
#'
#' \describe{
#' \item{pvalue}{The pvalue of the enrichment test.}
#' \item{lm.fit}{The linear model result. The pvalues are adjusted for ks.test}
#' \item{group}{The categrories of variables}
#' \item{hyper.matrix}{enrichment.hyper only. The hypergeometric test parameters
#' for each variable.}
#' \item{p.cutoff}{enrichment.hyper only. The pvalue cutoff to define
#' whether a variable was affected or not.}
#' \item{d}{enrichment.ks only. The value of ks.test statistics.
#' See \code{\link{ks.test}}}.
#' }
#' \item{adj.p.val}{enrichment.ks only. The adjusted pvalues for each variable}
#'
#' @export
#' @importFrom magrittr `%>%`
#'
#' @examples
#' data(glyco)
#' X = glyco$conc_table
#' design = model.matrix(~flipgroup, data = as(glyco$sample_table, "data.frame"))
#' group = glyco$feature_data$Protein
#' enrichment_test(X, design, group, "flipgroupcLNS", "hyper")
#' enrichment_test(X, design, group, "flipgroupcLNS", "ks")
#'
#' @author Chenghao Zhu
#' @seealso \code{\link{ks.test}} \code{\link{phyper}}
enrichment_test = function(X, design, group, coef, test, p.cutoff = 1){
    stopifnot(ncol(X) == nrow(design))
    stopifnot(nrow(X) == length(group))
    stopifnot(length(coef) == 1)
    stopifnot(coef %in% colnames(design))
    stopifnot(test %in% c("hyper", "ks"))
    stopifnot(between(p.cutoff, 0, 1))

    fits = lapply(seq_len(nrow(X)), function(i) {
        x = X[i,]
        fit = lm(x ~ design + 0)
        df = fit$df.residual
        coef = summary(fit)$coefficients[paste0("design",coef),]
        c(coef, df = df)
    }) %>%
        do.call(rbind, .) %>%
        `rownames<-`(rownames(X)) %>%
        as.data.frame %>%
        `colnames<-`(c("estimate", "stderr", "tvalue", "pvalue", "df"))

    if(test == "hyper"){
        res = .test_hyper(fits, group, p.cutoff)
    } else {
        res = .test_ks(fits, group)
    }
    return(res)
}

#' @importFrom magrittr `%>%`
#' @keywords internal
.test_hyper = function(fits, group, p.cutoff){
    res = lapply(c("greater", "less"), function(alt){
        compare = switch(
            alt,
            "greater" = function(x){x > 0},
            "less" = function(x){x < 0}
        )
        lapply(unique(group), function(ele){
            N = nrow(fits)
            m = sum(group == ele)
            n = N - m
            k = sum(compare(fits$estimate) & fits$pvalue < p.cutoff)
            x = sum(compare(fits$estimate) & fits$pvalue < p.cutoff & group == ele)
            p = phyper(q = x - 1, m = m, n = n, k = k, lower.tail = FALSE)
            return(c(N = N, m = m, n = n, k = k, x = x, p = p))
        }) %>%
            do.call(rbind, .) %>%
            `rownames<-`(unique(group))
    }) %>%
        `names<-`(c("greater", "less"))
    res = structure(
        list(
            pvalue = sapply(res, function(r) r[, "p"]),
            hyper.matrix = lapply(res, function(r) r[,1:5]),
            p.cutoff = p.cutoff,
            lm.fit = fits,
            group = group
        ),
        class = "enrichment.hyper"
    )
    return(res)
}

#' @importFrom magrittr `%>%`
#' @keywords internal
.test_ks = function(fits, group){
    res = lapply(c("greater", "less"), function(alt){
        pvalues = pt(fits$tvalue, fits$df, lower.tail = alt == "less")
        res = lapply(unique(group), function(ele){
            pvals = pvalues[group == ele]
            ks = ks.test(
                pvals, seq(from = 0, to = 1, length.out = length(pvals)),
                alternative = "greater"
            )
            d = ks$statistic
            names(d) = NULL
            c("d" = d, "p" = ks$p.value)
        }) %>%
            do.call(rbind,.) %>%
            `rownames<-`(unique(group))
        names(pvalues) = rownames(fits)
        list(pvalues = res[, "p"], "d" = res[, "d"], adjusted.p.values = pvalues)
    }) %>%
        `names<-`(c("greater", "less"))
    res = structure(
        list(
            pvalue = sapply(res, function(r) r$pvalues),
            d = sapply(res, function(r) r$d),
            lm.fit = fits,
            group = group,
            adj.p.val = sapply(res, function(r) r$adjusted.p.values)
        ),
        class = "enrichment.ks"
    )
    return(res)
}

#' @title enrichment barplot
#' @title This function makes bar plot for top most positively and negatively
#' enriched levels.
#' @param x enrichment.hyper. Returned by \code{\link{enrichement_test}}
#' when test was set to hyper.
#' @param each integer. The top most enriched levels on each side.
#' @return a ggplot object
#' @export
#' @import ggplot2
#' @import dplyr
#' @import cowplot
#' @author Chenghao Zhu
#' @seealso \code{\link{enrichment_test}}
plot.enrichment.hyper = function(x, each = 3, ...){
    df = x$hyper.matrix %>% as.data.frame
    df$var = rownames(df)
    df = df %>%
        mutate(diff.x = greater.x - less.x) %>%
        arrange(diff.x) %>%
        filter(var %in% var[1:each] | var %in% rev(var)[1:each]) %>%
        mutate(var = factor(var, levels = var))
    p1 = ggplot(df) +
        geom_col(aes(x = var, y = greater.x), fill = "gray10", width = .75) +
        geom_text(
            aes(x = var, y = greater.x + max(greater.x)/30, label = greater.x)
        ) +
        coord_flip() +
        labs(x = "") +
        theme_classic() +
        theme(axis.line.x = element_blank(),
              axis.text.x = element_blank(),
              axis.ticks.x = element_blank(),
              axis.title.x = element_blank(),
              axis.text.y = element_text(
                  color = "black", hjust = 0.5, margin = margin(0, 20, 0, 0)
              ),
              axis.title.y = element_blank())
    p2 = ggplot(df) +
        geom_col(aes(x = var, y = less.x), fill = "gray10", width = .75) +
        geom_text(aes(x = var, y = less.x + max(less.x)/30, label = less.x)) +
        scale_y_reverse() +
        scale_x_discrete(position = "top") +
        coord_flip() +
        labs(x = "") +
        theme_classic() +
        theme(axis.line.x = element_blank(),
              axis.text.x = element_blank(),
              axis.ticks.x = element_blank(),
              axis.title.x = element_blank(),
              axis.text.y = element_blank())
    plot_grid(p2, p1, ncol = 2, rel_widths = c(1, 1.1))
}

#' @title ecdf plot for enrichment test
#' @description This function makes a ecdf plot for enrichment test result
#' tested using ks.test.
#' @param x enrichment.ks. Returned by \code{\link{enrichement_test}}
#' when test was set to ks.
#' @param level character. The level to plot the ecdf plot. Must be in the
#' group.
#' @param alt character. Must be either greater or less.
#' @return a ggplot2 object
#' @export
#' @import ggplot2
#' @import dplyr
#' @author Chenghao Zhu
#' @seealso \code{\link{enrichment_test}}
plot.enrichment.ks = function(x, level, alt, ...){
    stopifnot(level %in% x$group)
    stopinfot(alt %in% colnames(x$pvalue))
    df = data.frame(
        pvalue = x$adj.p.val[,alt]
    )
    df$var = rownames(df)
    df %>%
        mutate(group = x$group) %>%
        filter(group == level) %>%
        ggplot() +
        stat_ecdf(geom = "step", aes(pvalue), color = "blue") +
        stat_ecdf(geom = "point", aes(pvalue), color = "blue") +
        stat_ecdf(geom = "step", aes(seq(0,1,length.out = length(pvalue))),
                  color = "red") +
        stat_ecdf(geom = "point", aes(seq(0,1,length.out = length(pvalue))),
                  color = "red") +
        labs(y = "Fn(pvalue)", title = level) +
        theme_classic() +
        theme(plot.title = element_text(hjust = 0.5))
}
