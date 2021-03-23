\name{modifyDataset}
\alias{modifyDataset}

\title{Modification of the values of the variables in the dataset}

\description{The function modifies the values of the declared variables used 
in the random part of the model if they are not unique. 
Unique values of the variables are required to 
build correct Z matrix for unsampled population elements.
}
\usage{
modifyDataset(data, names)
}

\arguments{
  \item{data}{the population dataset.
}
  \item{names}{the vector of names of the dataset columns which values should 
  be modified (names of the variables used to define the random part 
  of the model).
}
}
\value{
The dataset with modified values of the declared variables.
}

\author{
Tomasz Zadlo}

\examples{
data(realestData)
# some values of "NUTS2" and "NUTS4type" are the same - we will modify them:
modifyDataset(realestData, c("NUTS2", "NUTS4type")) 
} 


