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

output_consent <- consents |>
  filter(!is.na(fecha_ci)) |>
  transmute(
    id_paciente = nhc,
    ci = if_else(!is.na(fecha_ci), "Yes", "No"),
    ci_type = 1, # Adult
    ci_date = fecha_ci,
    ci_version = 6,
    ci_options__1 = 1,
    ci_options__2 = if_else(comparte_con_grupos, 1, NA),
    ci_options__3 = if_else(comparte_con_registros, 1, NA),
    ci_options__4 = if_else(consiente_relacionados, 1, NA),
    ci_opt_pers_data__0 = if_else(!consiente_contacto & is.na(email), 1, NA),
    ci_opt_pers_data__1 = if_else(!consiente_contacto & !is.na(email), 1, NA),
    ci_opt_pers_data__2 = if_else(consiente_contacto & !is.na(email), 1, NA),
    ci_options_assays__0 = if_else(!interesado_ensayos, 1, NA),
    ci_options_assays__1 = if_else(interesado_ensayos, 1, NA),
    ci_date_remove = fecha_retirada_ci,
    ci_update = pmax(fecha_ci, fecha_retirada_ci),
  )

output_personalinformation <- consents |>
  left_join(
    ufela_pacientes |>
      left_join(ufela_seguimiento, by = "pid") |>
      select(nhc, cip, fecha_ultima_visita, exitus, fecha_exitus),
    by = "nhc"
  ) |>
  filter(!is.na(fecha_ci) | exitus | !is.na(fecha_exitus)) |>
  transmute(
    id_paciente = nhc,
    pat_mail = email,
    nhc = cip,
    pat_referred = NA,
    hosp_refer_reason = NA,
    pat_derived = NA,
    hosp_deriv_reason = NA,
    pat_ref_hosp_follow = case_when(
      !is.na(fecha_exitus) ~ 0,
      !is.na(fecha_perdida_seguimiento) ~ 0,
      (fecha_ultima_visita - today()) <= dyears(1) ~ 1,
      TRUE ~ 98
    )
  )

  )
