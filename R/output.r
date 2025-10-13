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

output_demographics <- consents |>
  left_join(
    ufela_pacientes |> select(
      pid, nhc, sexo, fecha_nacimiento, provincia_nacimiento,
      codigo_postal, estudios, situacion_laboral_actual,
      situacion_laboral_actual_otra, exitus, fecha_exitus
    ),
    by = "nhc"
  ) |>
  left_join(ufela_clinica |> select(pid, fumador, starts_with("historia_familiar_")), by = "pid") |>
  transmute(
    family_history_dic = if_else(
      historia_familiar_alzheimer | historia_familiar_parkinson | historia_familiar_motoneurona,
      1, 0
    ),
    als_fam_hist_type = if_else(historia_familiar_motoneurona, 1, 2),
    adopted = 98,
    father_alive = NA,
    father_death_cause = NA,
    father_death_age = NA, # text?
    mother_alive = NA,
    mother_death_cause = NA,
    mother_death_age = NA, # text?
    fam_tree = 0,
    adv_fam = 0,
    sex = sexo |> case_match("Hombre" ~ 1, "Mujer" ~ 2, .default = 98),
    dob = fecha_nacimiento,
    cp_current = codigo_postal,
    res_type_current = NA,
    cp_life = NA,
    res_type_life = NA,
    country_origin = NA,
    province_origin = provincia_nacimiento,
    ehnicity_2_yn = NA,
    consanguinity = 98,
    consanguinity_type = NA,
    dominance = NA,
    smoke = case_match(fumador, 
      "No fumador" ~ 0,
      "Exfumador" ~ 1,
      "Fumador" ~ 2,
    ),
    start_smoke_age = NA,
    cigar_smoked = NA,
    end_smoke_age = NA,
    alcohol = NA,
    alcohol_freq = NA,
    alcohol_units = NA,
    start_alcohol_age = NA,
    end_alcohol_age = NA,
    adv_risk = 0,
    adv_med_history = 0,
    educ_level = estudios |> case_match(
      "No sabe leer ni escribir" ~ 0, # ISCED 0
      "Primarios incompletos" ~ 0, # ISCED 0
      "Primarios completos" ~ 1, # ISCED 1
      "Secundarios (ESO, BUP, COU, etc)" ~ 2, # ISCED 2
      "FP (grado medio o superior)" ~ 4, # ISCED 4-5
      "Universidad" ~ 6, # ISCED 6-7
      "Doctorado" ~ 8, # ISCED 8
      "Otro" ~ NA,
    ),
    educ_years = NA,
    employ_premorbid = NA,
    net_income_bef_onset = NA,
    income_periods = NA,
    employ_current = situacion_laboral_actual |> case_match(
      "Trabaja" ~ 1,
      "Parado" ~ 3,
      "Parado con subsidio | Prestación" ~ 3,
      "Labores de la casa" ~ 5,
      "Jubilado" ~ 6,
      "Incapacitado (o con invalidez permanente)" ~ 7,
      "Otra" ~ 99
    ),
    employ_current_other = situacion_laboral_actual_otra |> na_if(""),
    incapacity = coalesce(if_else(situacion_laboral_actual == 7, NA, 0), 98),
    incapacity_type = 98,
    disability = coalesce(if_else(situacion_laboral_actual == 7, NA, 0), 98),
    previous_jobs = 0,
  )
