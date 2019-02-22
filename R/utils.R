#' Make tibble (data frame) for one group
#'
#' Make tibble from one vector x. The vector is re-organised in a
#' data frame with column 1 for group label and column 2 for observations.
#'
#' @param x A vector.
#' @param gr_names The name of the first column - default "gr" for groups.
#' @param obs_names The name of the second column - default "obs" for
#'   observations.
#' @return A tibble - see details in \code{\link[TIBBLE]{tibble}}.
#' @examples
#' # Generate one vector
#' g1 <- rnorm(10)
#'
#' # Make tibble using default parameters
#' df <- mkt1(g1)
#'
#' # Make tibble using custom parameters
#' df <- mkt1(g1, gr_names = "groups", obs_names = "values", group_label =
#' "group1")
#'
#' # For instance, vector x:
#' x <- c(1, 2, 3, 4, 5)
#' df <- mkt1(x)
#' # become a tibble: 5 × 2:
#'        gr   obs
#'    <fctr> <dbl>
#' 1  Group1     1
#' 2  Group1     2
#' 3  Group1     3
#' 4  Group1     4
#' 5  Group1     5
#' @export
mkt1 <- function(x, name = c("obs", "gr")){
  gr <- c(rep('Group1', length(x)))
  df <- data.frame(x, gr)
  names(df) <- name
  df <- tibble::as_tibble(df)
  df
}

#' Make tibble (data frame) for two groups
#'
#' Make tibble from two vectors x and y. The two vectors are re-organised in a
#' data frame with column 1 for group labels as factors and column 2 for observations.
#'
#' @param x,y Two vectors.
#' @param gr_names The name of the first column - default "gr" for groups.
#' @param obs_names The name of the second column - default "obs" for
#'   observations.
#' @return A tibble - see details in \code{\link[TIBBLE]{tibble}}.
#' @examples
#' # Generate two vectors
#' g1 <- rnorm(10)
#' g2 <- rnorm(10)
#'
#' # Make tibble using default parameters
#' df <- mkt2(g1, g2)
#'
#' # Make tibble using custom parameters
#' df <- mkt2(g1, g2, gr_names = "groups", obs_names = "values", group_labels =
#' c("group1", "group2"))
#'
#' # For instance, vectors x & y:
#' x <- c(1, 2, 3, 4, 5)
#' y <- c(11, 12, 13, 14, 15)
#' df <- mkt2(x, y)
#' # become a tibble: 10 × 2:
#'        gr   obs
#'    <fctr> <dbl>
#' 1  Group1     1
#' 2  Group1     2
#' 3  Group1     3
#' 4  Group1     4
#' 5  Group1     5
#' 6  Group2    11
#' 7  Group2    12
#' 8  Group2    13
#' 9  Group2    14
#' 10 Group2    15
#' @export
mkt2 <- function(x, y, gr_names = "gr", obs_names = "obs",
                 group_labels = c("Group1", "Group2")){
  gr <- factor(c(rep(group_labels[1], length(x)),
          rep(group_labels[2], length(y)))) # group labels
  obs <- c(x, y) # observations
  df <- tibble::tibble(gr, obs) # make tibble
  names(df) <- c(gr_names, obs_names)
  df
}

subset_formula <- function(data, formula){
  vars <- all.vars(formula)
  param_col_name <- vars[2]
  obs_col_name <- vars[1]
  # check that the columns exist
  if (!(param_col_name %in% colnames(data))) {
    stop(paste0(param_col_name," does not exist"))
  }
  if (!(obs_col_name %in% colnames(data))) {
    stop(paste0(obs_col_name," does not exist"))
  }
  # check that param_col is a factor
  if(!is.factor(data[[param_col_name]])){
    stop('The predictor column must be a factor')
  }
  # check that obs_col is numeric
  if(!is.numeric(data[[obs_col_name]])){
    stop('Data must be numeric')
  }
  # get levels of param_col_name
  gr_names <- levels(data[[param_col_name]])
  # gr_name1 <- gr_names[[1]]
  # gr_name2 <- gr_names[[2]]
  # mm <- model.matrix(formula, data = data)
  # x <- data[mm[,2]==0,obs_col_name]
  # y <- data[mm[,2]==1,obs_col_name]
  # x <- x[[1]]
  # y <- y[[1]]
  # outputs
  out <- list(param_col_name = param_col_name,
              obs_col_name = obs_col_name,
              gr_names = gr_names)
}

subset_formula_wide <- function(data, formula){
  vars <- all.vars(formula)
  x_col_name <- vars[2]
  y_col_name <- vars[1]
  # check that the columns exist
  if (!(x_col_name %in% colnames(data))) {
    stop(paste0(x_col_name," does not exist"))
  }
  if (!(y_col_name %in% colnames(data))) {
    stop(paste0(y_col_name," does not exist"))
  }
  # check that param_col is numeric
  if(!is.numeric(data[[x_col_name]])){
    stop('The x column must be numeric')
  }
  # check that obs_col is numeric
  if(!is.numeric(data[[y_col_name]])){
    stop('The y column must be numeric')
  }
  # outputs
  out <- list(x_col_name = x_col_name,
    y_col_name = y_col_name)
}


elimna <- function(m){
  #
  # remove any rows of data having missing values
  #
  # From Rallfun-v32.txt - see \url{https://github.com/nicebread/WRS/}
  if(is.list(m)){
    for(j in 1:length(m))m[[j]]=na.omit(m[[j]])
    elimna=m
  }
  if(!is.list(m)){
    if(is.null(dim(m)))m<-as.matrix(m)
    ikeep<-c(1:nrow(m))
    for(i in 1:nrow(m))if(sum(is.na(m[i,])>=1))ikeep[i]<-0
    elimna<-m[ikeep[ikeep>=1],]
  }
  elimna
}
