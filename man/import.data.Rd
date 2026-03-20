\name{import.data}
\alias{importdata}
\alias{import_data}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Import .dta dataset
}
\description{
Imports .dta dataset; optionally, keeps columns you specify; optionally creates a date variables.
}
\usage{
import.data(file_path, ..., day = "none", month = "none", year = "none", name_date = "date")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{file_path}{
Path to file: see getfp() output
}
  \item{...}{
Set of columns you want kept in addition to the ID and, optionally, date column
}
  \item{day}{
Column name that specifies the "day" portion of the date; optional.
}
  \item{month}{
Column name that specifies the "day" portion of the date; optional.
}
  \item{year}{
Column name that specifies the "day" portion of the date; optional.
}
  \item{name_date}{
Name for the date column.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Returns a dataset
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
preg = import.data("PREGHIST.DTA", everything(), day = "daypregh", month = "monpregh", year = "yrpreghi", name_date = "dateintvw")
}
