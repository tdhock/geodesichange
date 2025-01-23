plot_model <- function(model_dt){
  if(requireNamespace("ggplot2")){
    ggplot2::ggplot()+
      ggplot2::theme_bw()+
      ggplot2::geom_vline(ggplot2::aes(
        xintercept=x),
        color="grey",
        data=result$model[
        , .SD[, .(x=unique(c(min_param,max_param)))]
        , by=data_i
        ])+
      ggplot2::geom_segment(ggplot2::aes(
        min_param, min_param*Linear+Constant,
        xend=max_param, yend=max_param*Linear+Constant),
        data=result$model)+
      ggplot2::facet_grid(data_i ~ .)+
      ggplot2::scale_x_continuous(breaks=seq(0,360,by=90))
  }
}
