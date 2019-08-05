#' Visualize Uncertainty in Projection
#'
#' This function simulates multivariate Poisson count data from a latent factor model.
#' @param coord1 - mcmc object with first coordinate of projection
#' @param coord2 - mcmc object with second coordinate of projection
#' @param type - plot type, options are 'points', 'circles', 'scatter'
#' @param items - items to display uncertainty, only works for circles or scatter
#' @return plot.obj - ggplot object
#' @importFrom magrittr "%>%"
#' @export

CredibleViz <- function(coord1, coord2, type = 'points', items = NULL){
  ### Create thin data
  num.pts <- ncol(coord1)

  coord1.thin <- coord1 %>% as.data.frame() %>% tidyr::gather %>% dplyr::mutate(key = ordered(key, levels=unique(key)))
  coord2.thin <- coord2 %>% as.data.frame() %>% tidyr::gather %>% dplyr::mutate(key = ordered(key, levels=unique(key)))
  combined <- data.frame(lv1 = coord1.thin$value, lv2 = coord2.thin$value, key = coord1.thin$key)
  combined.summarized <- combined %>% dplyr::group_by(key) %>% dplyr::summarise(mean.lv1 = mean(lv1), mean.lv2 = mean(lv2)) %>% dplyr::mutate(id = 1:num.pts)

  ### Points plot
  plot.obj <- ggplot2::ggplot(combined.summarized, aes(x=mean.lv1, y=mean.lv2)) + ggplot2::geom_text(aes(label=id, color=as.factor(id))) + ggplot2::guides(color=FALSE) + ggplot2::xlab('latent variable one') + ggplot2::ylab('latent variable two')  + ggplot2::theme_bw() + ggplot2::ggtitle('')

  gg_color_hue <- function(n) {
    hues = seq(15, 375, length = n + 1)
    hcl(h = hues, l = 65, c = 100)[1:n]
  }
  cols = gg_color_hue(num.pts)

  if (type == 'scatter'){
    ### Scatter plot
    # match colors

    for (item.numb in 1:length(items)){
      sims.tmp <- combined %>% dplyr::filter(key == items[item.numb])
      item.pos <- (1:num.pts)[colnames(coord1) == items[item.numb]]
      plot.obj <- plot.obj + ggplot2::geom_point(data = sims.tmp, aes(x=lv1, y=lv2), alpha=.15, color = cols[item.pos])
      plot.obj
    }
    plot.obj <- plot.obj + ggplot2::geom_text(aes(label=id, color=as.factor(id))) + ggplot2::ylim(-3.1,3.1) + ggplot2::xlim(-3.1,3.1)
    plot.obj
  } else if (type == 'circles'){

    ### circle plot
    all.itempos <- (1:num.pts)[colnames(coord1) %in% items]

    for (item.numb in 1:length(items)){
      sims.temp <- combined %>% dplyr::filter(key == items[item.numb]) %>% dplyr::select(lv1,lv2) %>% as.matrix
      fit <- MASS::cov.mve(sims.temp, quantile.used = nrow(sims.temp) * 0.95, nsamp = 100)
      points_in_ellipse <- sims.temp[fit$best, ]
      ellipse_boundary <- predict(cluster::ellipsoidhull(points_in_ellipse))
      ellipse_boundary_df <- data.frame(ellipse_boundary)
      plot.obj <- plot.obj + ggplot2::geom_path(data=ellipse_boundary_df, aes(x=V1, y=y),linetype=3, color=cols[all.itempos[item.numb]])
    }

    plot.obj <- plot.obj + ggplot2::ylim(-3.1,3.1) + ggplot2::xlim(-3.1,3.1)
  }

  return(list(plot.obj = plot.obj))
}

