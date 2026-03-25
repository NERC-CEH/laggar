#' Assembling
#'
#' @param response Vector of response values
#' @param covariate Matrix of covariate values, assembled by lg_ summary functions
#' @param index_seq Sequence of index values used to index the covariate matrix
#' @param minindex Minimum index, ignored if \code{index_seq} specified
#' @param maxindex Maximum index, ignored if \code{index_seq} specified
#' @param incindex index increment, ignored if \code{index_seq} specified
#'
#' @returns list
#' @export
#'
lg_assembly <- function(response, covariate,
                        index_seq = NULL,
                        minindex = NULL, maxindex = NULL, incindex = NULL){
  index <- .index_check(minindex = minindex, maxindex = maxindex,
                        incindex = incindex, index_seq = index_seq)
  if(length(index) != ncol(covariate)){
    warning(paste("Index sequence must be of the same length as the number of",
                  "columns within the covariate matrix. Replacing the supplied",
                  "index with a sequence from 1 to",ncol(covariate)))
    index <- seq(1,ncol(covariate))
  }
  # checks
  # check if response is a numeric vector
  if(!class(response) %in% c("integer","numeric"))
    stop("response must be a numeric or integer vector")
  # check if covariate is a matrix
  if(!("matrix" %in% class(covariate)))
    stop("covariate must be a matrix")

  if(!identical(nrow(covariate), length(response)))
    stop("number of rows of covariate must be equal to length of response")
  # check for NAs
  if(anyNA(covariate)){
    cov_na_index <- apply(covariate, 1, anyNA)
    if(all(cov_na_index)){
      stop("covariate matrix has NAs present in every row")
    } else{
      message("Removing rows where covariate matrix has NAs")
      cov_notna_index <- !cov_na_index
    }
  } else{
    cov_notna_index <- rep(TRUE, length(response))
  }
  if(anyNA(response)){
    resp_na_index <- is.na(response)
    if(all(resp_na_index)){
      stop("Response is entirely NAs")
    } else{
      message("Removing rows where response is NA")
      resp_notna_index <- !resp_na_index
    }
  }else{
    resp_notna_index <- rep(TRUE, length(response))
  }
  notna_index <- (cov_notna_index+resp_notna_index)==2
  response <- response[notna_index]
  covariate <- covariate[notna_index,]

  index_matrix <- matrix(index,
                        nrow = nrow(covariate),
                        ncol = ncol(covariate),
                        byrow = TRUE)

  data_list <- list(response = response,
                    covariate = covariate,
                    index = index_matrix)

  return(data_list)

}



### Helper functions ########################

.sf_conversion <- function(object, name = "object"){
  if(!(any(c("sf","sfc","sfg") %in% class(object)))){
    message(paste("Converting",name,"to sf object"))
    object <- sf::st_as_sf(object)
  }
  return(object)
}


.index_check <- function(minindex,maxindex,incindex,index_seq){
  if(any(c(is.null(minindex),is.null(maxindex),is.null(incindex))) & is.null(index_seq))
    stop("Either index_seq or minindex, maxindex and incindex must be specified")
  if(is.null(index_seq)){
    if(!all(is.numeric(minindex),is.numeric(maxindex),is.numeric(incindex)))
      stop("minindex, maxindex and incindex must be numeric")
    if(!all(length(minindex)==1,length(maxindex)==1,length(incindex)==1))
      stop("minindex, maxindex and incindex must be of length 1")
    if(minindex <= 0){
      stop("minindex must be greater than 0")
    }
    if(maxindex < minindex){
      stop("maxindex must be greater than minindex")
    }
    index <- seq(minindex,maxindex,incindex)
  } else{
    if(!is.numeric(index_seq))
      stop("index_seq must be numeric")
    if(length(index_seq) < 2)
      stop("index_seq must be of length greater than 1")
    if(min(index_seq)<= 0){
      stop("Smallest indexance increment must be greater than 0")
    }
    index <- index_seq
  }
  return(index)
}
