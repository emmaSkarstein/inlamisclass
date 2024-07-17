#' Birthweight data
#'
#' A dataset with 189 observations, with the birth weight of a child as the response, and the smoking status of the mother, as well as the mother's weight at the time of the last menstrual cycle as covariates.
#'
#' @format ## `birthweight`
#' A data frame with 189 rows and 3 columns:
#' \describe{
#'   \item{bwt}{birth weight of the child (in kg, centered)}
#'   \item{lwt}{mother's weight at last menstrual cycle (in kg, centered)}
#'   \item{smoke}{mother's smoking status (0 = not smoking, 1 = smoking), suspected to be misreported}
#' }
#' @source Hosmer et al. (2013), <https://books.google.no/books?id=64JYAwAAQBAJ>
"birthweight"

#' Case-control study on cervical cancer
#'
#' A dataset with 2044 observations, with whether or not a patient has cervical cancer as response, and exposure to herpes simplex virus type 2 measured by two different procedures as covariates. The first procedure is the most accurate, and that measurement is recorded as `x`, whereas the other procedure is less accurate and is recorded as `w`. The more accurate measurement is only available for 115 patients, while the less accurate measurement is available for all subjects, meaning that we have a validation data set.
#'
#' @format ## `case_control_data`
#' A data frame with 189 rows and 3 columns:
#' \describe{
#'   \item{d}{cervical cancer status (1 = cases (has cervical cancer), 0 = controls)}
#'   \item{x}{measurement from a refined western blot procedure (more accurate)}
#'   \item{w}{measurement from less accurate western blot procedure}
#' }
#' @source Carroll, Gain and Lubin (1993), <https://doi.org/10.1080/01621459.1993.10594310>
"case_control_data"
