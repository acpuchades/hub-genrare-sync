library(dplyr)
library(stringr)
library(lubridate)

source(here::here("R", "consents.r"))
source(here::here("R", "precision.r"))
source(here::here("R", "ufela.r"))

parse_phenotype_information <- function(phenotype, phenotype_other, als=1, umn=2, lmn=3, fosmn=4) {
  case_match(phenotype,
      "ELA Espinal" ~ als,
      "ELA Bulbar" ~ als,
      "ELA Respiratoria" ~ als,
      "Monomiélica" ~ als,
      "Pseudopolineurítica" ~ als,
      "Hemipléjica (Mills)" ~ als,
      "Parálisis bulbar progresiva" ~ als,
      "Esclerosis Lateral Primaria (ELP)" ~ umn,
      "Atrofia Muscular Progresiva (AMP)" ~ lmn,
      "Otro" ~ case_match(phenotype_other |> str_to_upper(),
        "ELA ESPINAL-DLFT" ~ als,
        "ELA ESPINAL-DLFT-PK" ~ als,
        "ELA BULBAR-DLFT" ~ als,
        "ELA BULBAR-DEMENCIA" ~ als,
        "ELA ESPINAL-PARKINSONISMO" ~ als,
        "ELA BULBAR +/- DLFT" ~ als,
        "ELA ESPINAL POSIBLE" ~ als,
        "ELA ESPINAL VS HIPOPARATIROIDISMO" ~ als,
        "ELA  EXPANSION C9ORF72" ~ als,
        "ELA FAMILIAR" ~ als,
        "ELA GENERALIZADA" ~ als,
        "ELA INICIO RESPIRATORIO" ~ als,
        "ENM 1+2" ~ als,
        "ELADFT" ~ als,
        "ELA+PARKINSONISMO" ~ als,
        "ELA-DFT" ~ als,
        "ESPINAL-DLFT-PK" ~ als,
        "PSEUDOPOLIRRADICULAR" ~ als,
        "MN Y PK" ~ als,
        "PBP" ~ als,
        "ELP VS HSP" ~ umn,
        "MALALTIA SEGONA MOTONEURONA" ~ lmn,
        "FORMA DIFUSA DE SEGUNDA MOTONEURONA ASOCIADA A ANTI-GM1" ~ lmn,
        "SEGUNDA MOTONEURONA DIFUSA" ~ lmn,
        "FOSMN" ~ fosmn,
      )
    )
}

redcap_site_id <- "1234"

output_patients_ci <- consents |>
  filter(!is.na(fecha_ci)) |>
  select(nhc) |>
  drop_na()

output_patients_dead <- ufela_pacientes |>
  filter(exitus | !is.na(fecha_exitus)) |>
  left_join(ufela_clinica |> select(pid, starts_with("fenotipo_")), by = "pid") |>
  filter(fenotipo_al_diagnostico != "Otro" | !str_detect(fenotipo_al_diagnostico_otro,
    "AME|ATAXIA|ATROFIA MUSCULAR ESPINAL|FEWDON|IBP|KENNEDY|MG|NO ELA|PARAPARESIA|POSTPOLIO"
  )) |>
  select(nhc) |>
  drop_na()

output_patient_ids_allocated <- tibble()

output_patient_ids_unallocated <- output_patients_ci |>
  full_join(output_patients_dead, by = "nhc") |>
  mutate(record_id = str_glue("{redcap_site_id}-{row_number()}"))

output_patient_ids <- output_patient_ids_allocated |>
  bind_rows(output_patient_ids_unallocated) |>
  left_join(ufela_pacientes |> select(nhc, pid, cip), by = "nhc")

output_status <- output_patient_ids |>
  left_join(consents, by = "nhc") |>
  inner_join(ufela_pacientes |> select(nhc, exitus, fecha_exitus), by = "nhc") |>
  transmute(
    record_id,
    ci_date_entry = fecha_ci,
    status = case_when(
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

output_consent <- output_patient_ids |>
  inner_join(consents, by = "nhc") |>
  transmute(
    record_id,
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

output_personalinformation <- output_patient_ids |>
  left_join(
    consents |> select(nhc, email, fecha_perdida_seguimiento),
    by = "nhc"
  ) |>
  left_join(
    ufela_pacientes |>
      left_join(ufela_seguimiento, by = "pid") |>
      select(nhc, fecha_ultima_visita, exitus, fecha_exitus),
    by = "nhc"
  ) |>
  transmute(
    record_id,
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

output_demographics <- output_patient_ids |>
  left_join(
    ufela_pacientes |> select(
      pid, sexo, fecha_nacimiento, provincia_nacimiento,
      codigo_postal, estudios, situacion_laboral_actual,
      situacion_laboral_actual_otra, exitus, fecha_exitus
    ),
    by = "pid"
  ) |>
  left_join(ufela_clinica |> select(pid, fumador, starts_with("historia_familiar_")), by = "pid") |>
  transmute(
    record_id,
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

output_clinicaldata <- output_patient_ids |>
  left_join(
    ufela_clinica |> select(
      pid, fecha_inicio_clinica, fecha_diagnostico_ELA,
      fenotipo_al_diagnostico, fenotipo_al_diagnostico_otro
    ),
    by = "pid"
  ) |>
  left_join(
    ufela_nutri |> select(
      pid, fecha_basal_nutri = "fecha_visita", estatura, peso, imc_actual, peso_premorbido, fecha_peso_premorbido
    ),
    by = join_by(pid, closest(fecha_diagnostico_ELA <= fecha_basal_nutri))
  ) |>
  left_join(
    suppressWarnings(ufela_alsfrs |> summarize(
      alsfrs_data_available = any(!is.na(total)),
      kings_data_available = any(!is.na(kings)),
      fecha_perdida_lenguaje = pick(fecha_visita, lenguaje) |>
        filter(lenguaje == 0) |>
        pull(fecha_visita) |>
        min(na.rm = TRUE) |>
        na_if(as_date(Inf)),
      fecha_perdida_deambulacion = pick(fecha_visita, caminar) |>
        filter(caminar <= 1) |>
        pull(fecha_visita) |>
        min(na.rm = TRUE) |>
        na_if(as_date(Inf)),
      fecha_vmni_22h = pick(fecha_visita, insuficiencia_respiratoria) |>
        filter(insuficiencia_respiratoria <= 1) |>
        pull(fecha_visita) |> 
        min(na.rm = TRUE) |>
        na_if(as_date(Inf)),
      fecha_iot = pick(fecha_visita, insuficiencia_respiratoria) |>
        drop_na() |>
        filter(
          insuficiencia_respiratoria == 0,
          last(insuficiencia_respiratoria == 0)
        ) |>
        pull(fecha_visita) |>
        min(na.rm = TRUE) |>
        na_if(as_date(Inf)),
      .by = pid
    )),
    by = "pid"
  ) |>
  left_join(
    suppressWarnings(ufela_respi |> summarize(
      fecha_indicacion_vmni = min(fecha_indicacion_vmni, na.rm = TRUE) |> na_if(as_date(Inf)),
      .by = pid
    )),
    by = "pid"
  ) |>
  left_join(
    suppressWarnings(ufela_nutri |> summarize(
      fecha_indicacion_gastrostomia = min(fecha_indicacion_peg, na.rm = TRUE) |> na_if(as_date(Inf)),
      .by = "pid"
    )),
    by = "pid"
  ) |>
  left_join(
    pals_ecas |> summarize(
      ecas_data_available = any(!is.na(fecha)),
      .by = "nhc"
    ),
    by = "nhc"
  ) |>
  left_join(
    pals_eq5 |> summarize(
      eq5_data_available = any(!is.na(fecha)),
      .by = "nhc"
    ),
    by = "nhc"
  ) |>
  transmute(
    record_id,
    date_clinical_onset = fecha_inicio_clinica,
    simmetry = NA,
    laterality_onset = NA,
    als_symp_bulbar__3 = if_else(fenotipo_al_diagnostico %in% c("ELA Bulbar", "Parálisis bulbar progresiva"), 1, 0),
    als_symp_spinal__6 = if_else(fenotipo_al_diagnostico == "Flail arm", 1, 0),
    als_symp_spinal__7 = if_else(fenotipo_al_diagnostico == "Flail leg", 1, 0),
    als_symp_resp__1 = if_else(fenotipo_al_diagnostico == "ELA respiratoria", 1, 0),
    # orpha_code,
    als_symp_other_old = fenotipo_al_diagnostico,
    height = estatura,
    weight_at_diagnosis = peso,
    weight_reference = peso_premorbido,
    weight_reference_date = fecha_peso_premorbido,
    metabolic = 0,
    vital_signs = 0,
    clinical_phenotype = parse_phenotype_information(fenotipo_al_diagnostico, fenotipo_al_diagnostico_otro, fosmn = NA),
    strength = 0,
    neuromus_data = 0,
    speech_loss = if_else(!is.na(fecha_perdida_lenguaje), 1, 0),
    speech_loss_date = fecha_perdida_lenguaje,
    ambulation_loss = if_else(!is.na(fecha_perdida_deambulacion), 1, 0),
    ambulation_loss_date = fecha_perdida_deambulacion,
    tracheostomy_ind = if_else(!is.na(fecha_iot), 1, 0),
    tracheostomy_ind_date = fecha_iot,
    niv_indicated_date = fecha_indicacion_vmni,
    niv_indicated = if_else(!is.na(fecha_indicacion_vmni), 1, 0),
    niv_22h = if_else(!is.na(fecha_vmni_22h), 1, 0),
    niv_22h_date = fecha_vmni_22h,
    gast_carrier_ind = if_else(!is.na(fecha_indicacion_gastrostomia), 1, 0),
    gast_ind_date = fecha_indicacion_gastrostomia,
    alsfrs_data = if_else(alsfrs_data_available, 1, 0) |> replace_na(0),
    cognitive_data = if_else(ecas_data_available, 1, 0) |> replace_na(0),
    kings_data = if_else(kings_data_available, 1, 0) |> replace_na(0),
    eq_5d_5l_data = if_else(eq5_data_available, 1, 0) |> replace_na(0),
    cgic_data = 0,
    pgic_data = 0,
  )

output_diagnosis <- output_patient_ids |>
  left_join(
    ufela_pacientes |> select(pid, exitus, fecha_exitus),
    by = "pid"
  ) |>
  left_join(
    ufela_clinica |> select(
      pid, fecha_inicio_clinica, fecha_diagnostico_ELA,
      fenotipo_al_diagnostico, fenotipo_al_diagnostico_otro,
      fenotipo_al_exitus, fenotipo_al_exitus_otro,
      starts_with("resultado_estudio_"),
    ),
    by = "pid"
  ) |>
  left_join(
    ufela_respi |> summarize(
      fvc_data_available = any(!is.na(fvc_sentado)),
      adv_resp_data_available = TRUE,
      .by = pid,
    ),
    by = "pid"
  ) |>
  transmute(
    record_id,
    initial_diagnosis = parse_phenotype_information(fenotipo_al_diagnostico, fenotipo_al_diagnostico_otro),
    date_diagnosis = fecha_diagnostico_ELA,
    final_diagnosis = case_when(
      initial_diagnosis %in% c(1, 4) ~ 1,
      initial_diagnosis %in% 2:3 ~ if_else(
        exitus & (fecha_exitus - fecha_inicio_clinica) < dyears(4),
        1, parse_phenotype_information(fenotipo_al_exitus, fenotipo_al_exitus_otro)
      )
    ),
    eer_category = NA,
    gold_coast_crit_yn = if_else(initial_diagnosis %in% c(1, 3), 1, NA),
    fvc_dic = if_else(fvc_data_available, 1, 0) |> replace_na(0),
    respiratory_data = if_else(adv_resp_data_available, 1, 0) |> replace_na(0),
    nfl_data = 98,
    other_liquid_biomarkers = 98,
    blood_analysis = 98,
    emg = 98,
    tms = 0,
    pet = 0,
    rm = 98,
    genetic_data = !is.na(resultado_estudio_c9) | !is.na(resultado_estudio_sod1),
    samples_data = 98,
  )

output_treatment <- output_patient_ids |>
  left_join(ufela_clinica |> select(pid, riluzol), by = "pid") |>
  left_join(
    ufela_respi |> summarize(
      uso_niv = any(portador_vmni, na.rm = TRUE),
      .by = pid,
    ),
    by = "pid"
  ) |>
  left_join(
    ufela_alsfrs |> 
      select(pid, fecha_visita, insuficiencia_respiratoria) |>
      drop_na() |>
      filter(insuficiencia_respiratoria == 0) |>
      summarize(
        fecha_iot = min(fecha_visita, na.rm = TRUE) |> na_if(as_date(Inf)),
        .by = pid
      ),
    by = "pid"
  ) |>
  left_join(
    ufela_nutri |>
      filter(portador_peg) |>
      slice_min(fecha_visita, by = pid, n = 1, na_rm = TRUE) |>
      select(pid, fecha_gastrostomia = "fecha_colocacion_peg"),
    by = "pid"
  ) |>
  left_join(
    ufela_nutri |> summarize(
      uso_supl_nutricional = any(suplementacion_nutricional_oral | suplementacion_nutricional_entera, na.rm = TRUE),
      .by = pid,
    ),
    by = "pid"
  ) |>
  transmute(
    record_id,
    treatment_dic = if_else(riluzol, 1, 0) |> replace_na(98),
    other_treatment_dic = 98,
    clinical_trials_dic = 98,
    niv_dic = if_else(uso_niv, 1, 0) |> replace_na(0),
    cough = 98,
    tracheostomy = if_else(!is.na(fecha_iot), 1, 0) |> replace_na(0),
    tracheostomy_date = fecha_iot,
    gast_carrier_dic = if_else(!is.na(fecha_gastrostomia), 1, 0) |> replace_na(0),
    gast_carrier_type = if_else(!is.na(fecha_gastrostomia), 0, NA), # PEG
    gast_date = fecha_gastrostomia,
    nut_sup = if_else(uso_supl_nutricional, 1, 0) |> replace_na(98),
  )

output_alstreatmentdata <- output_patient_ids |>
  left_join(ufela_clinica |> select(pid, riluzol, fecha_inicio_riluzol), by = "pid") |>
  filter(riluzol) |>
  transmute(
    record_id,
    treat_type = 1, # riluzole
    med_start_date = fecha_inicio_riluzol,
    treat_status = 1,
  )

output_alsfrs <- output_patient_ids |>
  inner_join(
    ufela_alsfrs |>
      select(pid, fecha_visita, lenguaje:insuficiencia_respiratoria, total, mitos) |>
      filter(if_any(c(lenguaje:insuficiencia_respiratoria, total, mitos), ~!is.na(.))),
    by = "pid", relationship="one-to-many"
  ) |>
  left_join(
    ufela_nutri |> select(pid, fecha_visita_nutri = "fecha_visita", fecha_colocacion_peg),
    by = join_by(pid, closest(fecha_visita >= fecha_visita_nutri)), multiple = "first"
  ) |>
  transmute(
    record_id,
    alsfrs_date = fecha_visita,
    alsfrs_total_man = total,
    mitos_staging_man = mitos,
    alsfrs_gast_use = case_when(
      is.na(fecha_colocacion_peg) ~ 0,
      fecha_visita < fecha_colocacion_peg ~ 0,
      fecha_visita >= fecha_colocacion_peg ~ 1,
    ),
    alsfrs_language = lenguaje,
    alsfrs_salivation = salivacion,
    alsfrs_deglution = deglucion,
    alsfrs_writing = escritura,
    alsfrs_cutting_nopeg = if_else(alsfrs_gast_use == 0, cortar_sin_peg, NA),
    alsfrs_cutting_peg__4 = if_else(alsfrs_gast_use == 1, cortar_con_peg, NA),
    alsfrs_clothing = vestido,
    alsfrs_bed = cama,
    alsfrs_walking = caminar,
    alsfrs_stairs = subir_escaleras,
    alsfrs_dyspnea = disnea,
    alsfrs_orthopnea = ortopnea,
    alsfrs_resp_failure = insuficiencia_respiratoria,
  )

output_kings <- output_patient_ids |>
  inner_join(
    ufela_alsfrs |> select(pid, fecha_visita, kings),
    by = "pid", relationship = "one-to-many"
  ) |>
  drop_na(fecha_visita, kings) |>
  transmute(
    record_id, 
    kings_date = fecha_visita,
    kings_total = kings
  )

output_eq5d5l <- output_patient_ids |>
  inner_join(pals_eq5, by = "nhc", relationship = "one-to-many") |>
  transmute(
    record_id,
    eq_5d_5l_date = fecha,
    mobility = movilidad,
    selfcare = autocuidado,
    usual_act = actividades,
    pain = dolor_malestar,
    anxiety = ansiedad_depresion,
    health_status = eva
  )
