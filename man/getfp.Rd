\name{getfp}
\alias{getfilepath}
\alias{file_path}
\alias{filepath}
\title{Get file path}
\usage{
getfp(file_name = "", memade = FALSE)
}
\description{
Combines the path of folder where you store your files and name of file to return path to file.
Path should be set in the "data_folder" variable in the costants.R file. If you are using the function to output a file to the output folder, set memade = TRUE.
}
\examples{
getfp(file_name = "blood.dta")
# returns "~/data/blood.dta"
}
