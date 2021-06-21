#' @import parsnip
#' @importFrom "tibble" "as_tibble"
#'
register_glasp_cox <- function(){

  try(set_new_model("glasp_cox"), silent=TRUE)

  set_model_engine("glasp_cox", mode = "classification", eng = "glasp")
  set_dependency("glasp_cox", eng = "glasp", pkg = "glasp")

  set_model_arg(
    model = "glasp_cox",
    eng = "glasp",
    parsnip = "l1",
    original = "l1",
    func = list(pkg = "dials", fun = "penalty", range = c(-6, 1)),
    has_submodel = F
  )

  set_model_arg(
    model = "glasp_cox",
    eng = "glasp",
    parsnip = "l2",
    original = "l2",
    func = list(pkg = "dials", fun = "penalty", range = c(-6, 1)),
    has_submodel = F
  )

  set_model_arg(
    model = "glasp_cox",
    eng = "glasp",
    parsnip = "frob",
    original = "frob",
    func = list(pkg = "dials", fun = "penalty", range = c(-6, 1)),
    has_submodel = F
  )

  set_model_arg(
    model = "glasp_cox",
    eng = "glasp",
    parsnip = "num_comp",
    original = "num_comp",
    func = list(pkg = "dials", fun = "num_comp", range = c(1L, 20L)),
    has_submodel = F
  )

  set_fit(
    model = "glasp_cox",
    eng = "glasp",
    mode = "classification",
    value = list(
      interface = "formula",
      protect = c("formula", "data"),
      func = c(pkg = "glasp", fun = "cox_regression"),
      defaults = list()
    )
  )

  prob_info <-  pred_value_template(
    post = function(x, object) {
      as_tibble(x)
    },
    func = c(fun = "predict"),
    object = quote(object$fit),
    new_data = quote(new_data),
    type = "prob"
  )
  set_pred(
    model = "glasp_cox",
    eng = "glasp",
    mode = "classification",
    type = "prob",
    value = prob_info
  )

  num_info <-
    list(
      pre = NULL,
      post = NULL,
      func = c(fun = "predict"),
      args =
        # These lists should be of the form:
        # {predict.mda argument name} = {values provided from parsnip objects}
        list(
          # We don't want the first two arguments evaluated right now
          # since they don't exist yet. `type` is a simple object that
          # doesn't need to have its evaluation deferred.
          object = quote(object$fit),
          new_data = quote(new_data),
          type = "numeric"
        )
    )
  set_pred(
    model = "glasp_cox",
    eng = "glasp",
    mode = "classification",
    type = "numeric",
    value = num_info
  )

  class_info <-
    list(
      pre = NULL,
      post = NULL,
      func = c(fun = "predict"),
      args =
        # These lists should be of the form:
        # {predict.mda argument name} = {values provided from parsnip objects}
        list(
          # We don't want the first two arguments evaluated right now
          # since they don't exist yet. `type` is a simple object that
          # doesn't need to have its evaluation deferred.
          object = quote(object$fit),
          new_data = quote(new_data),
          type = "class"
        )
    )
  set_pred(
    model = "glasp_cox",
    eng = "glasp",
    mode = "classification",
    type = "class",
    value = num_info
  )

  set_encoding(
    model = "glasp_cox",
    eng = "glasp",
    mode = "classification",
    options = list(
      predictor_indicators = "none",
      compute_intercept = FALSE,
      remove_intercept = FALSE,
      allow_sparse_x = FALSE
    )
  )
}


#' parsnip interface for GLASP - Cox Models
#'
#' @param mode 	A single character string for the type of model.
#' @param l1 Regularization in norm l1 (lasso)
#' @param l2 Regularization in norm l2 (ridge)
#' @param frob Regularization in Frobenius norm. It is unique to glasp and controls the importance
#' of the clustering on the variables. If frob=0, we want to fit a model equivalent to elastic-net,
#' without clustering. If frob > l1+l2 means that we are more interested on the feature clustering
#' than the variable selection. If frob > 1 we are more interested on finding clusters than on
#' getting a linear model with good predictive power.
#' @param num_comp Maximum number of clusters to search.
#'
#' @return A `glasp_cox` parsnip model
#'
#' @export
#'
glasp_cox <- function(mode = "classification", l1=NULL, l2=NULL, frob=NULL, num_comp=NULL){
  args <- list(l1 = rlang::enquo(l1),
               l2 = rlang::enquo(l2),
               frob = rlang::enquo(frob),
               num_comp = rlang::enquo(num_comp))

  # Save some empty slots for future parts of the specification
  out <- list(args = args, eng_args = NULL,
              mode = mode, method = NULL, engine = NULL)

  # set classes in the correct order
  class(out) <- make_classes("glasp_cox")
  out
}
