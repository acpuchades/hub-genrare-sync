library(dplyr)
library(lubridate)

source(here::here("R", "consents.r"))
source(here::here("R", "ufela.r"))

output_status <- consents |>
  inner_join(ufela_pacientes |> select(nhc, cip, exitus, fecha_exitus), by = "nhc") |>
  filter(!is.na(fecha_ci) | exitus | !is.na(fecha_exitus)) |>
  transmute(
    id_paciente = nhc,
    ci_date_entry = fecha_ci,
    status = dplyr::case_when(
      exitus | !is.na(fecha_exitus) ~ 0, # Deceased
      !is.na(fecha_perdida_seguimiento) ~ 2, # Lost follow-up
      !is.na(fecha_retirada_ci) ~ 3, # Withdraws participation
      TRUE ~ 1 # Alive
    ),
    lost_date = fecha_perdida_seguimiento,
    ci_optedout_date = fecha_retirada_ci,
    death_date = fecha_exitus,
    status_related_exitus = NA,
    death_place = if_else(exitus, 98, NA_integer_), # Unknown
    status_update = if_else(exitus, fecha_exitus, today()),
    seed_project = if_else(!is.na(seed_als), "Yes", "No")
  )

  )
