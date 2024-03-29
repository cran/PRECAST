\name{plot_scatter}
\alias{plot_scatter}
\title{Scatter plot for two-dimensional embeddings}
\description{
  Scatter plot for two-dimensional embeddings
}
\usage{
  plot_scatter(embed_use, meta_data, label_name, 
    xy_names=c('tSNE1', 'tSNE2'), no_guides = FALSE, 
    cols = NULL, 
    point_size = 0.5, point_alpha=1, 
    base_size = 12, do_points = TRUE, do_density = FALSE, border_col='gray',
    legend_pos='right', legend_dir='vertical', nrow.legend=NULL)
}
\arguments{
  \item{embed_use}{an object named "Seurat", "maxtrix" or "dgCMatrix". The object of class "Seurat" must include slot "scale.data".}
  \item{meta_data}{an optional positive integer, specify the number of features to be extracted.}
  \item{label_name}{the size of point in the scatter plot.}
  \item{xy_names}{the text size in the plot.}
  \item{no_guides}{whether display the legend.}
  \item{cols}{colors used in the plot.}
  \item{point_size}{the point size of scatter plot.}
  \item{point_alpha}{the transparency of the plot.}
  \item{base_size}{the base text size.}
  \item{do_points}{Plot point plot.}
  \item{do_density}{Plot density plot}
  \item{border_col}{the border color in the plot.}
  \item{legend_pos}{the position of legend.}
  \item{legend_dir}{the direction of legend.}
  \item{nrow.legend}{the number of rows of legend.}
}
\details{
  Nothing
}
\value{
  Return a ggplot2 object.
}

\author{
  Wei Liu
}
\note{
  nothing
}
\seealso{
  None
}
\examples{
  embed_use <- cbind(1:100, sin((1:100)*pi/2))
  meta_data <- data.frame(cluster=factor(rep(1:2, each=50)))
  plot_scatter(embed_use, meta_data, label_name='cluster')
}