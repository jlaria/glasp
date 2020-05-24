#' @import parsnip
#' @importFrom "tibble" "as_tibble"
#'
.onLoad <- function(libname, pkgname){
  try(set_new_model("glasp_model"), silent = TRUE)


  # Linear regression -------------------------------------------------------
  set_model_engine(
    "glasp_model",
    mode = "regression",
    eng = "glasp"
  )
  set_dependency("glasp_model", eng = "glasp", pkg = "glasp")

  set_model_arg(
    model = "glasp_model",
    eng = "glasp",
    parsnip = "l1",
    original = "l1",
    func = list(pkg = "dials", fun = "penalty", range = c(-6, 1)),
    has_submodel = F
  )
  set_model_arg(
    model = "glasp_model",
    eng = "glasp",
    parsnip = "l2",
    original = "l2",
    func = list(pkg = "dials", fun = "penalty", range = c(-6, 1)),
    has_submodel = F
  )
  set_model_arg(
    model = "glasp_model",
    eng = "glasp",
    parsnip = "frob",
    original = "frob",
    func = list(pkg = "dials", fun = "penalty", range = c(-6, 1)),
    has_submodel = F
  )
  set_model_arg(
    model = "glasp_model",
    eng = "glasp",
    parsnip = "num_comp",
    original = "num_comp",
    func = list(pkg = "dials", fun = "num_comp", range = c(1L, 20L)),
    has_submodel = F
  )

  set_fit(
    model = "glasp_model",
    eng = "glasp",
    mode = "regression",
    value = list(
      interface = "formula",
      protect = c("formula", "data"),
      func = c(pkg = "glasp", fun = "linear_regression"),
      defaults = list()
    )
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
    model = "glasp_model",
    eng = "glasp",
    mode = "regression",
    type = "numeric",
    value = num_info
  )


  # Classification ----------------------------------------------------------
  set_model_engine(
    "glasp_model",
    mode = "classification",
    eng = "glasp"
  )

  set_fit(
    model = "glasp_model",
    eng = "glasp",
    mode = "classification",
    value = list(
      interface = "formula",
      protect = c("formula", "data"),
      func = c(pkg = "glasp", fun = "linear_classification"),
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
    model = "glasp_model",
    eng = "glasp",
    mode = "classification",
    type = "prob",
    value = prob_info
  )

}


# glasp_model -------------------------------------------------------------


#' parsnip interface for Variable Cluster - Linear Regression Models
#'
#' @param mode 	A single character string for the type of model.
#' @param l1 An non-negative number representing the total amount of regularization L1
#' @param l2 An non-negative number representing the total amount of group regularization L2
#' @param frob An non-negative number representing the total amount of regularization in the Frobenius norm of the error of the low-rank approximation
#' @param num_comp Number of clusters to search for
#'
#' @export
#' @examples
#' \dontrun{
#' library(glasp)
#' set.seed(0)
#' # generate data -----------------------------------------------------------
#' data <- glasp::simulate_dummy_surv_data(N = 100)
#' # compute glasp -----------------------------------------------------------
#' glasp_mod <- glasp::glasp_model(l1 = tune::tune(),
#'                                 l2 = tune::tune(),
#'                                 frob = tune::tune(),
#'                                 num_comp = 3)
#'
#' glasp_mod <- parsnip::set_mode(glasp_mod, "classification")
#' glasp_mod <- parsnip::set_engine(glasp_mod, "glasp")
#'
#' data_rs <- rsample::vfold_cv(data, v = 10)
#' metric_vals <- yardstick::metric_set(yardstick::roc_auc)
#' ctrl <- tune::control_grid(verbose = T)
#'
#' rec_form <- tune::tune_grid(glasp_mod, event~.,
#'                             resamples = data_rs,
#'                             metrics = metric_vals,
#'                             control = ctrl,
#'                             grid = 100)
#'
#' best_config <- show_best(rec_form, metric = "yardstick::roc_auc")
#' model <- glasp::linear_classification( event~., data,
#'                                        l1 = best_config$l1[1],
#'                                        l2 = best_config$l2[1],
#'                                        frob = best_config$frob[1],
#'                                        num_comp = 3)
#'
#' model$beta
#' }
glasp_model <-
  function(mode = "unknown",  l1 = NULL, l2 = NULL, frob = NULL, num_comp = NULL) {
    # Check for correct mode

    # Capture the arguments in quosures
    args <- list(l1 = rlang::enquo(l1),
                 l2 = rlang::enquo(l2),
                 frob = rlang::enquo(frob),
                 num_comp = rlang::enquo(num_comp))

    # Save some empty slots for future parts of the specification
    out <- list(args = args, eng_args = NULL,
                mode = mode, method = NULL, engine = NULL)

    # set classes in the correct order
    class(out) <- make_classes("glasp_model")
    out
  }
