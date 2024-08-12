#' Plotting for DLIM-IM
#' @description Create cumulative and point-wise effect plots for a DLIM-IM
#' @seealso \link[dlimIM]{pred_m}
#' @seealso \link[dlimIM]{MCMC_sampler_m}
#' @export
#' @import ggplot2
#' @param x posterior samples  (class "\code{dlimIM}")
#' @param pred prediction object from \link[dlimIM]{pred_m}
#' @param m_star vector of modifier index values (class "\code{numeric}")
#' @param burnin number of MCMC samples to remove as warm-up (class "\code{numeric}")
#' @param thin number post-burn-in MCMC samples to thin by (class "\code{numeric}")
#' @param type plot type options: "cumulative", "by_time", "by_modifier" (class "\code{character}")
#' @param exp_times labels for exposure time points (class "\code{character}")
#' @param time_pts exposure time points to plot by when \code{plot_by = "by_time"} (class "\code{numeric}")
#' @param n_col number of columns for plotting grid (class "\code{numeric}")
#' @return This function returns a ggplot of specified \code{type}


plot_bdlim <- function(x, pred = NULL, m_star, burnin,
                       thin = 1, type,
                       exp_times = NULL, time_pts=NULL,
                       n_col = 3){

  if(is.null(pred)){
    pred <- pred_m(x,
                   burnin = burnin,
                   thin = thin,
                   m_star = m_star)
  }

  if(type == "cumulative"){
    plot_df <- data.frame(Modifiers = c(m_star),
                          Cumul_Effect = c(pred$betas_cumul),
                          LB = c(pred$cumul_LB),
                          UB = c(pred$cumul_UB))

    ggplot(plot_df, aes(x=Modifiers,y=Cumul_Effect)) +
      geom_hline(yintercept = 0) +
      geom_ribbon(aes(ymin=LB, ymax=UB), alpha=0.5 , fill = "grey70")+
      geom_line()+
      xlab("Weighted Modifier Index") +
      ylab("Cumulative Effect") +
      theme_classic()
  }else if(type == "by_time" | type == "by_modifier"){
    if(type == "by_time"){
      X_lab <- "Weighted Modifier Index"
      xaxis <- m_star #put modifiers on x-axis
      plot_by <- time_pts #make plots by specific time points
      effect <- c(pred$betas[, which(exp_times %in% time_pts)])
      lb <- c(pred$LB[, which(exp_times %in% time_pts)])
      ub <- c(pred$UB[, which(exp_times %in% time_pts)])
    }else{
      X_lab <- "Exposure Time"
      xaxis <- exp_times
      plot_by <- m_star
      effect <- c(t(pred$betas))
      lb <- c(t(pred$LB))
      ub <- c(t(pred$UB))
    }

    plot_df <- data.frame(Xaxis = rep(xaxis, length(plot_by)),
                          Plot_by = factor(rep(plot_by, each = length(xaxis))),
                          Effect = effect,
                          LB = lb,
                          UB = ub)

    ggplot(plot_df, aes(x=Xaxis,y=Effect)) +
      geom_hline(yintercept = 0)+
      geom_ribbon(aes(ymin=LB, ymax=UB), alpha=0.5, color=FALSE)+
      geom_line()+
      facet_wrap(vars(Plot_by), ncol = n_col) +
      xlab(X_lab) +
      ylab("Linear effect per exposure unit") +
      theme_classic()

  }


}
