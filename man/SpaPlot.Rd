\name{SpaPlot}
\alias{SpaPlot}
\title{Spatial heatmap}
\description{
  Plot spatial heatmap for a Seurat object with spatial transcriptomics data.
}
\usage{
  SpaPlot(seuInt, batch=NULL, item=NULL, point_size=2,text_size=12, 
                    cols=NULL,font_family='', border_col="gray10",
                    fill_col='white', ncol=2, combine = TRUE,
                    title_name="Sample", ...)
}
\arguments{
  \item{seuInt}{an object named "Seurat".}
  \item{batch}{an optional positive integer or integer vector, specify the batches to be extracted. Users can check the batches' names by \code{unique(seuInt$batch)}.}
  \item{item}{an optional string, which column is plotted in the meta data of seuInt. Users can check the meta data by \code{head(seuInt@meta.data)}. If \code{item} takes value from ("RGB_UMAP", "RGB_tSNE"), this function will plot the RGB plot.}
  \item{point_size}{the size of point in the scatter plot.}
  \item{text_size}{the text size in the plot.}
  \item{cols}{colors used in the plot}
  \item{font_family}{the font family used for the plot, default as Times New Roman.}
  \item{border_col}{the border color in the plot.}
  \item{fill_col}{the  color used in backgroup.}
  \item{ncol}{the number of columns in the layout of plots.}
  \item{combine}{an optional logical value, whether plot all on a figure. If TRUE, all figures are plotted; otherwise, return a list with each plot as component.}
  \item{title_name}{an optional string, title name in the plot.}
  \item{...}{other arguments passed to \code{\link{plot_scatter}}}.
}
\details{
  Nothing
}
\value{
  Return a ggplot2 object or list of ggplots objects.
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

  data(PRECASTObj)
  PRECASTObj <- SelectModel(PRECASTObj)
  seuInt <- IntegrateSpaData(PRECASTObj, species='unknown')
  SpaPlot(seuInt)
  
}