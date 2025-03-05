#' Low Back Pain Diagnoses and Diagnosis Criteria
#' 
#' In a 2011 study, Smart et al. (2011) collected information on 464 Irish and
#' British patients suffering from low back pain regarding:
#' - the type of low back pain (classified into 'nociceptive', 'peripheral
#'   neuropathic', and 'central neuropathic')
#' - the presence/absence of 38 clinical criteria and symptoms relating to low
#'   back pain.
#'   
#' Fop et al. (2017) conducted Latent Class Analysis on this data set to
#' retrieve the experts' classifications; and by comparing nested models they
#' were able to select 11 out of 38 criteria which contain the most of the
#' relevant grouping information while avoiding redundancy.
#' 
#' @usage data('lowbackpain')
#' 
#' @format
#' A list containing:
#' \describe{
#'  \item{data:}{A 464x11 binary matrix indicating the presence/absence of the
#'              11 selected criteria for each of the 464 patients. Rownames are
#'              study IDs.}
#'  \item{group:}{A vector of length 464 indicating the diagnosis each patient
#'           received, numerically coded (order has no meaning). The names of
#'           the vector give the diagnosis in words.}
#'  \item{index:}{The index to the criteria explaining which symptom they refer to.}
#' }
#' 
#' @references
#' - Fop, M, Smart, K, Murphy, TB (2017).
#'   *Variable Selection for Latent Class Analysis with Application to Low Back Pain Diagnosis*.
#'   The Annals of Applied Statitics. 11(4), 2080-2110.
#'   \doi{10.1214/17-aoas1061}
#' - Smart, K, Blake, C, Staines, A, Doody, C (2011).
#'   *The Discriminative Validity of "Nociceptive", "Peripheral Neuropathic", and "Central Sensitization" as Mechanisms-Based Classifications of Musculoskeletal Pain*.
#'   The Clinical Journal of Pain. 27, 655-663.
#'   \doi{10.1097/AJP.0b013e318215f16a}
#'   
#' @source
#'
#' Supplemental Content for Fop et al. (2017): \doi{10.1214/17-AOAS1061SUPP}
"lowbackpain"
