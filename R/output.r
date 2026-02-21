library(dplyr)
library(here)
library(lubridate)
library(purrr)
library(readr)
library(stringi)
library(stringr)
library(tidyr)

source("R/biobanco.r")
source("R/consents.r")
source("R/precision.r")
source("R/ufela.r")

#redcap_site_id <- "1196"
redcap_site_id <- "1939"

output_tofersen_study_ids <- c(
  "1196-6", "1196-7", "1196-8", "1196-9"
)

data_dir <- here("data", "20260221")
genrare_patient_ids_path <- file.path(data_dir, "genrare-patient-ids.csv")

get_redcap_site_id <- function(x) {
  str_extract(x, "^([0-9]+)-([0-9]+)$", group = 1) |> as.integer()
}

get_redcap_patient_id <- function(x) {
  str_extract(x, "^([0-9]+)-([0-9]+)$", group = 2) |> as.integer()
}

make_redcap_table <- function(data, event_name,
                              record_id = "record_id",
                              repeat_instrument = NULL)
{
  res <- data |>
    mutate(across(where(is.timepoint), ~strftime(.x, "%Y-%m-%d"))) |>
    mutate(redcap_event_name = event_name, .after = all_of(record_id))

  if (is.null(repeat_instrument)) {
    res |>
      arrange(
        get_redcap_site_id(.data[[record_id]]),
        get_redcap_patient_id(.data[[record_id]])
      )
  } else {
    res |>
      mutate(
        redcap_repeat_instrument = repeat_instrument,
        .after = redcap_event_name
      ) |>
      mutate(
        redcap_repeat_instance = row_number(),
        .by = all_of(record_id), .after = redcap_repeat_instrument
      ) |>
      arrange(
        get_redcap_site_id(.data[[record_id]]),
        get_redcap_patient_id(.data[[record_id]]),
        redcap_repeat_instance
      )
  }
}

write_redcap_table <- function(data, file) {
  write_csv(data, file, na = "")
}

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

output_patient_ids_allocated <- readr::read_csv(genrare_patient_ids_path) |>
  janitor::clean_names()

output_patients_ci <- consents |>
  drop_na(nhc, fecha_ci) |>
  select(nhc)

output_patients_dead <- ufela_pacientes |>
  filter(exitus | !is.na(fecha_exitus)) |>
  left_join(ufela_clinica |> select(pid, starts_with("fenotipo_")), by = "pid") |>
  filter(fenotipo_al_diagnostico != "Otro" | !str_detect(fenotipo_al_diagnostico_otro,
    "AME|ATAXIA|ATROFIA MUSCULAR ESPINAL|FEWDON|KENNEDY|MG|NO ELA|PARAPARESIA|POSTPOLIO"
  )) |>
  select(nhc) |>
  drop_na()

redcap_last_record_id <- coalesce(
  output_patient_ids_allocated |>
    pull(record_id) |>
    str_extract("-([0-9]+)$", group=1) |>
    as.integer() |>
    max(), 0
)

output_patient_ids_unallocated <- output_patients_dead |>
  full_join(output_patients_ci, by = "nhc") |>
  anti_join(output_patient_ids_allocated, by = "nhc") |>
  left_join(
    ufela_pacientes |>
      select(nhc, cip, pid, created_datetime),
    by = "nhc"
  ) |>
  arrange(created_datetime) |>
  mutate(record_id = str_glue("{redcap_site_id}-{redcap_last_record_id+row_number()}")) |>
  select(record_id, nhc, cip, pid)

output_patient_ids <- output_patient_ids_allocated |>
  bind_rows(output_patient_ids_unallocated) |>
  arrange(
    get_redcap_site_id(record_id),
    get_redcap_patient_id(record_id)
  )

output_status <- output_patient_ids |>
  left_join(consents, by = "nhc") |>
  left_join(
    ufela_pacientes |>
      select(nhc, exitus, fecha_nacimiento, fecha_exitus),
    by = "nhc"
  ) |>
  transmute(
    record_id,
    ci_date_entry = fecha_ci,
    status = case_when(
      exitus | !is.na(fecha_exitus) ~ 0, # Deceased
      !is.na(fecha_perdida_seguimiento) ~ 2, # Lost follow-up
      !is.na(fecha_retirada_ci) ~ 3, # Withdraws participation
      TRUE ~ 1 # Alive
    ),
    dob = fecha_nacimiento,
    lost_date = fecha_perdida_seguimiento,
    ci_optedout_date = fecha_retirada_ci,
    death_date = fecha_exitus,
    status_related_exitus = NA,
    death_place = if_else(exitus, 98, NA_integer_), # Unknown
    status_update = today(),
    seed_project = if_else(!is.na(seed_als), 1, 0),
    status_complete = if_else(
      is.na(exitus) | exitus,
      1, # Unverified (missing death info)
      2  # Complete
    )
  )

output_consent <- output_patient_ids |>
  inner_join(
    consents |> filter(!is.na(fecha_ci)),
    by = "nhc"
  ) |>
  transmute(
    record_id,
    ci = 1,
    ci_type = 1, # Adult
    ci_date = fecha_ci,
    ci_version = 6,
    ci_options___1 = 1,
    ci_options___2 = if_else(comparte_con_grupos, 1, 0),
    ci_options___3 = if_else(comparte_con_registros, 1, 0),
    ci_options___4 = if_else(consiente_relacionados, 1, 0),
    ci_opt_pers_data___0 = if_else(!consiente_contacto &  is.na(email), 1, 0),
    ci_opt_pers_data___1 = if_else(!consiente_contacto & !is.na(email), 1, 0),
    ci_opt_pers_data___2 = if_else( consiente_contacto & !is.na(email), 1, 0),
    ci_options_assays___0 = if_else(!interesado_ensayos, 1, 0),
    ci_options_assays___1 = if_else( interesado_ensayos, 1, 0),
    ci_date_remove = fecha_retirada_ci,
    ci_update = today(),
    consent_complete = 2, # Complete
  ) |>
  arrange(record_id |> str_extract("^([0-9]+)-([0-9]+)$", group=2) |> as.integer())

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
    pat_ref_hosp_follow = NA,
    personal_update = today(),
    personal_information_complete = 2, # Complete
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
  left_join(
    ufela_clinica |>
      select(pid, fumador, starts_with("historia_familiar_")),
    by = "pid"
  ) |>
  left_join(
    pals_patients |>
      select(nhc, pals_visita_inclusion = "visita_inclusion"),
    by = "nhc"
  ) |>
  left_join(
    pals_ecas |>
      select(nhc, ecas_education_years = "escolaridad"),
    by = "nhc"
  ) |>
  transmute(
    record_id,
    family_history_dic = if_else(historia_familiar_motoneurona, 1, 0),
    als_fam_hist_type = if_else(historia_familiar_motoneurona, 1, 2),
    adopted = 98,
    father_alive = NA,
    father_death_cause = NA,
    father_death_age = NA,
    mother_alive = NA,
    mother_death_cause = NA,
    mother_death_age = NA,
    fam_tree = 0,
    adv_fam = if_else(!is.na(pals_visita_inclusion), 1, 0),
    sex = sexo |> case_match("Hombre" ~ 1, "Mujer" ~ 2, .default = 98),
    cp_current = if_else(
      codigo_postal |> str_detect("^[0-9]+$"), # assumes ES type zipcodes
      codigo_postal, NA_character_
    ),
    res_type_current = NA,
    cp_life = NA,
    res_type_life = NA,
    country_origin = NA,
    province_origin = provincia_nacimiento,
    ethnicity = NA,
    ethnicity_2_yn = NA,
    ethnicity_2 = NA,
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
    adv_med_history = if_else(!is.na(pals_visita_inclusion), 1, 0),
    marital = NA,
    convivence_status = NA,
    educ_level = coalesce(
      case_when(
        is.na(ecas_education_years) ~ NA,
        ecas_education_years < 6 ~ 0,  # < Primaria
        ecas_education_years < 10 ~ 1, # Primaria
        ecas_education_years < 13 ~ 3, # ESO / FP medio
        ecas_education_years < 15 ~ 5, # FP superior
        ecas_education_years < 17 ~ 6, # Grado
        ecas_education_years < 19 ~ 7, # Máster
        ecas_education_years >= 19 ~ 8 # Doctorado
      ),
      estudios |> case_match(
        "No sabe leer ni escribir" ~ 0, # ISCED 0
        "Primarios incompletos" ~ 0, # ISCED 0
        "Primarios completos" ~ 1, # ISCED 1
        "Secundarios (ESO, BUP, COU, etc)" ~ 3, # ISCED 3
        "FP (grado medio o superior)" ~ 3, # ISCED 3/5
        "Universidad" ~ 6, # ISCED 6-7
        "Doctorado" ~ 8, # ISCED 8
        "Otro" ~ NA,
      )
    ),
    educ_old = NA,
    educ_years = ecas_education_years,
    employ_premorbid = NA,
    net_income_bef_onset = NA,
    income_periods = NA,
    employ_current = situacion_laboral_actual |> case_match(
      "Trabaja" ~ 1,
      "Parado" ~ 3,
      "Parado con subsidio | Prestación" ~ 3,
      "Estudiante" ~ 4,
      "Labores de la casa" ~ 5,
      "Jubilado" ~ 6,
      "Incapacitado (o con invalidez permanente)" ~ 7,
      "Otra" ~ 99
    ),
    employ_current_other = NA,
    incapacity = if_else(situacion_laboral_actual == 7, 1, 0),
    incapacity_type = 98,
    perm_incapacity_date = NA,
    disability = if_else(situacion_laboral_actual == 7, 1, 0),
    disability_perc = NA,
    dependency = 98,
    dependency_grade = 98,
    previous_jobs = 0,
    demo_update = today(),
    demographics_complete = case_when(
      !is.na(pals_visita_inclusion) ~ 1, # Unverified (missing pals data)
      situacion_laboral_actual_otra != "" ~ 1, # Unverified (check working status)
      historia_familiar_motoneurona ~ 1, # Unverified (check fALS status)
      TRUE ~ 2 # Complete
    ),
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
      pid, fecha_basal_nutri = "fecha_visita", estatura, peso,
      imc_actual, peso_premorbido, fecha_peso_premorbido
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
    als_symp_bulbar___3 = if_else(fenotipo_al_diagnostico %in% c("ELA Bulbar", "Parálisis bulbar progresiva"), 1, 0),
    als_symp_spinal___6 = if_else(fenotipo_al_diagnostico == "Flail arm", 1, 0),
    als_symp_spinal___7 = if_else(fenotipo_al_diagnostico == "Flail leg", 1, 0),
    als_symp_spinal___8 = if_else(fenotipo_al_diagnostico %in% c(
      "ELA Espinal", "Monomiélica", "Hemipléjica (Mills)", "Pseudopolineurítica"
    ), 1, 0),
    als_symp_resp___3 = if_else(fenotipo_al_diagnostico == "ELA Respiratoria", 1, 0),
    orpha_code = NA,
    als_symp_other_old___1 = if_else(fenotipo_al_diagnostico == "Otro", 1, 0),
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
    clinical_update = today(),
    clinical_data_complete = 2, # Complete
  )

output_geneticstudy <- bind_rows(
  output_patient_ids |>
    inner_join(ufela_clinica |> select(pid, starts_with("estudio_genetico_")), by = "pid") |>
    filter(estudio_genetico_c9 | estudio_genetico_atxn2 | estudio_genetico_ar) |>
    transmute(
      record_id,
      analy_seq = 6, # Expansion (triple PCR)
      analy_gene_triplepcr___ar = if_else(estudio_genetico_ar, 1, 0),
      analy_gene_triplepcr___atxn2 = if_else(estudio_genetico_atxn2, 1, 0),
      analy_gene_triplepcr___c9orf72 = if_else(estudio_genetico_c9, 1, 0),
    ),
  output_patient_ids |>
    inner_join(ufela_clinica |> select(pid, starts_with("estudio_genetico_")), by = "pid") |>
    filter(estudio_genetico_sod1 | estudio_genetico_fus | estudio_genetico_tardbp | estudio_genetico_unc13a) |>
    transmute(
      record_id,
      analy_seq = 1, # Sanger or Reduced panel
      analy_gene_sanger___sod1 = if_else(estudio_genetico_sod1, 1, 0),
      analy_gene_sanger___tardbp = if_else(estudio_genetico_tardbp, 1, 0),
      analy_gene_sanger___fus = if_else(estudio_genetico_fus, 1, 0),
      analy_gene_sanger___unc13a = if_else(estudio_genetico_unc13a, 1, 0),
    )
) |>
  mutate(
    genetstudy_update = today(),
    genetic_study_complete = 2 # Complete
  )

output_geneticvariations <- bind_rows(
  output_patient_ids |>
    filter(!(record_id %in% output_tofersen_study_ids)) |>
    inner_join(ufela_clinica |> select(pid, resultado_estudio_c9), by = "pid") |>
    filter(resultado_estudio_c9 == "Alterado") |>
    transmute(
      record_id,
      gen_gene = "HGNC:28337", # C9orf72
      gen_alle1_var_typ_2 = 2, # Repetitions
      gen_alle1_cnv_n = 999, # Repeated expansion
      gen_alle1_fx = 4, # Pathogenic
      clin_criteria = 1, # Yes
      genetic_variations_complete = 2, # Complete
    ),
  output_patient_ids |>
    filter(!(record_id %in% output_tofersen_study_ids)) |>
    inner_join(
      ufela_clinica |>
        select(pid, resultado_estudio_atxn2, starts_with("numero_repeticiones_atxn2_")),
      by = "pid"
    ) |>
    filter(resultado_estudio_atxn2 == "Alterado") |>
    transmute(
      record_id,
      gen_gene = "HGNC:10555", # ATXN2
      gen_alle1_var_typ_2 = 2, # Repetitions
      gen_alle1_cnv_n = numero_repeticiones_atxn2_a1,
      gen_alle2_cnv_n = numero_repeticiones_atxn2_a2,
      clin_criteria = 1, # Yes
      genetic_variations_complete = 2, # Complete
    ),
  output_patient_ids |>
    filter(!(record_id %in% output_tofersen_study_ids)) |>
    inner_join(ufela_clinica |> select(pid, resultado_estudio_sod1), by = "pid") |>
    filter(resultado_estudio_sod1 == "Alterado") |>
    transmute(
      record_id,
      gen_gene = "HGNC:11179", # SOD1
      gen_alle1_var_typ_2 = 1, # SNV
      clin_criteria = 98, # Unknown
      genetic_variations_complete = 0, # Incomplete
    ),
  output_patient_ids |>
    filter(!(record_id %in% output_tofersen_study_ids)) |>
    inner_join(ufela_clinica |> select(pid, resultado_estudio_tardbp), by = "pid") |>
    filter(resultado_estudio_tardbp == "Alterado") |>
    transmute(
      record_id,
      gen_gene = "HGNC:11571", # TARDBP
      gen_alle1_var_typ_2 = 1, # SNV
      clin_criteria = 98, # Unknown
      genetic_variations_complete = 0, # Incomplete
    ),
  output_patient_ids |>
    filter(!(record_id %in% output_tofersen_study_ids)) |>
    inner_join(ufela_clinica |> select(pid, resultado_estudio_fus), by = "pid") |>
    filter(resultado_estudio_fus == "ALTERADO") |>
    transmute(
      record_id,
      gen_gene = "HGNC:4010", # FUS
      gen_alle1_var_typ_2 = 1, # SNV
      clin_criteria = 98, # Unknown
      genetic_variations_complete = 0, # Incomplete
    ),
  output_patient_ids |>
    filter(!(record_id %in% output_tofersen_study_ids)) |>
    inner_join(ufela_clinica |> select(pid, resultado_estudio_unc13a), by = "pid") |>
    filter(resultado_estudio_unc13a == "Alterado") |>
    transmute(
      record_id,
      gen_gene = "HGNC:23150", # UNC13A
      gen_alle1_var_typ_2 = 1, # SNV
      clin_criteria = 98, # Unknown
      genetic_variations_complete = 0, # Incomplete
    ),
) |>
  mutate(genetanaly_update = today())

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
  left_join(
    bind_rows(
      biobanco_muestras_ela |> select(nhc),
      biobanco_muestras_pals |> select(nhc)
    ) |>
      summarize(
        samples_data_available = TRUE,
        .by = "nhc"
      ),
    by = "nhc"
  ) |>
  transmute(
    record_id,
    initial_diagnosis = parse_phenotype_information(fenotipo_al_diagnostico, fenotipo_al_diagnostico_otro),
    date_diagnosis = fecha_diagnostico_ELA,
    final_diagnosis = case_when(
      initial_diagnosis == 1 ~ 1,
      initial_diagnosis %in% 2:3 ~ if_else(
        exitus & (fecha_exitus - fecha_inicio_clinica) < dyears(4),
        1, parse_phenotype_information(fenotipo_al_exitus, fenotipo_al_exitus_otro)
      ),
      initial_diagnosis == 4 ~ NA,
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
    genetic_data = if_else(!is.na(resultado_estudio_c9) | !is.na(resultado_estudio_sod1), 1, 0) |> replace_na(0),
    samples_data = if_else(samples_data_available, 1, 0) |> replace_na(0),
    diagnosis_update = today(),
    diagnosis_complete = case_when(
      initial_diagnosis == 1 ~ 2, # Complete
      initial_diagnosis %in% c(2, 3) ~ if_else(
        (fecha_exitus - fecha_inicio_clinica) < dyears(4),
        2, # Complete (final_diagnosis=als)
        1  # Unverified (check final phenotype)
      ),
      initial_diagnosis == 4 ~ 1, # Unverified (check final diagnosis)
      TRUE ~ 2, # Complete
    )
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
    ufela_nutri |>
      summarize(
        uso_supl_nutricional = any(
          suplementacion_nutricional_oral |
            suplementacion_nutricional_entera,
          na.rm = TRUE
        ),
        .by = pid
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
    trmt_update = today(),
    treatment_complete = 2 # Complete
  )

output_alstreatmentdata <- output_patient_ids |>
  left_join(
    ufela_clinica |> select(pid, riluzol, fecha_inicio_riluzol),
    by = "pid"
  ) |>
  filter(riluzol) |>
  transmute(
    record_id,
    treat_type = 1, # riluzole
    med_start_date = fecha_inicio_riluzol,
    treat_status = 1,
    als_treatment_data_complete = 2 # Complete
  ) |>
  mutate(medication_update = today())

output_alsfrs <- output_patient_ids |>
  inner_join(
    ufela_alsfrs |>
      select(pid, fecha_visita, lenguaje:insuficiencia_respiratoria, total, mitos) |>
      filter(if_any(c(lenguaje:insuficiencia_respiratoria, total, mitos), ~!is.na(.))),
    by = "pid", relationship="one-to-many"
  ) |>
  left_join(
    ufela_nutri |>
      select(pid, fecha_visita_nutri = "fecha_visita", fecha_colocacion_peg),
    by = join_by(pid, closest(fecha_visita >= fecha_visita_nutri)),
    multiple = "first"
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
    alsfrs_cutting_peg = if_else(alsfrs_gast_use == 1, cortar_con_peg, NA),
    alsfrs_clothing = vestido,
    alsfrs_bed = cama,
    alsfrs_walking = caminar,
    alsfrs_stairs = subir_escaleras,
    alsfrs_dyspnea = disnea,
    alsfrs_orthopnea = ortopnea,
    alsfrs_resp_failure = insuficiencia_respiratoria,
    alsfrs_update = today(),
    alsfrs_complete = 2 # Complete
  )

output_kingsscale <- output_patient_ids |>
  inner_join(
    ufela_alsfrs |> select(pid, fecha_visita, kings),
    by = "pid", relationship = "one-to-many"
  ) |>
  drop_na(fecha_visita, kings) |>
  transmute(
    record_id,
    kings_date = fecha_visita,
    kings_total = kings,
    kings_update = today(),
    kings_scale_complete = 2 # Complete
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
    health_status = eva,
    hrqol_update = today(),
    eq_5d_5l_complete = 2, # Complete
  )

output_fvcinfo <- output_patient_ids |>
  inner_join(
    ufela_respi |>
      select(pid, fecha_visita, starts_with("fvc_")) |>
      filter(if_any(starts_with("fvc_"), ~!is.na(.))),
    by = "pid", relationship = "one-to-many"
  ) |>
  transmute(
    record_id,
    fvc_date = fecha_visita,
    fvc_basal_perc = fvc_sentado,
    fvc_basal_ml = fvc_sentado_absoluto,
    fvc_decub_perc = fvc_estirado,
    fvc_decub_ml = fvc_estirado_absoluto,
    fvc_update = today(),
    fvc_info_complete = 2 # Complete
  )

output_advancedrespiratorydata <- output_patient_ids |>
  inner_join(
    ufela_respi |>
      select(
        pid, fecha_visita,
        pim, pem, pns, # Respiratory strength
        ph_sangre_arterial, paco2, pao2, hco3, # Blood gases
        sao2_media, ct90, # Night pulsioximetry
      ) |>
      filter(if_any(-c(pid, fecha_visita), ~!is.na(.))),
    by = "pid", relationship = "one-to-many"
  ) |>
  transmute(
    record_id,
    respiratory_date = fecha_visita,
    mip_cmh2o = pim,
    mep_cmh2o = pem,
    snip = pns,
    sleep_under90 = ct90,
    spo2_media = sao2_media,
    ph = ph_sangre_arterial,
    pco2 = paco2,
    po2 = pao2,
    hco3 = hco3,
    resp_update = today(),
    advanced_respiratory_data_complete = 2 # Complete
  )

output_weightandbmi <- output_patient_ids |>
  inner_join(
    ufela_nutri |>
      select(pid, fecha_visita, peso) |>
      drop_na(),
    by = "pid", relationship = "one-to-many"
  ) |>
  transmute(
    record_id,
    weight_current_date = fecha_visita,
    weight_current = peso,
    bmi_update = today(),
    weight_and_bmi_complete = 2 # Complete
  )

output_cognitiveassessment <- output_patient_ids |>
  inner_join(pals_ecas, by = "nhc", relationship = "one-to-many") |>
  full_join(
    pals_alsftdq |> rename(alsftdq_total = total),
    by = c("nhc", "pals_id", "fecha")
  ) |>
  drop_na(record_id, fecha) |>
  transmute(
    record_id,
    cognitive_date = fecha,

    # ECAS
    nps_cognitive_test___1 = if_else(version != "-", 1, 0),
    ecas_version = case_match(version, "A" ~ 1, "B" ~ 2, "C" ~ 3),
    ecas_language = idioma |> case_match("Español" ~ 24),
    ecas_modality = modo |> case_match("Oral" ~ 1, "Escrito" ~ 2),
    ecas_language_naming = nombrar,
    ecas_language_comprehension = comprension,
    ecas_memory_immediate = recuerdo,
    ecas_language_spelling = deletreo,
    ecas_fluency_letter_p = fluidez_p,
    ecas_verbal_fluency_index = NA,
    ecas_executive_reverse_digits = secuencia,
    ecas_executive_alternation = alternancia,
    ecas_fluency_letter_t = fluidez_t,
    ecas_visuospatial_dots = puntos,
    ecas_visuospatial_cubes = cubos,
    ecas_visuospatial_numbers = numeros,
    ecas_executive_sentences = frases,
    ecas_social_cognition = social,
    ecas_memory_delayed_recall = retencion,
    ecas_memory_delayed_recognition = reconocimiento,
    ecas_language_total_man = lenguaje,
    ecas_verbal_fluency_total_man = fluidez,
    ecas_executive_total_man = ejecutiva,
    ecas_memory_total_man = memoria,
    ecas_visuospatial_total_man = visuoespacial,
    ecas_als_specific_man = ela_e,
    ecas_als_non_specific_man = ela_ne,
    ecas_total_man = total,

    # Behavioural (ECAS-CI)
    nps_cognitive_test___3 = if_else(!is.na(total_ec), 1, 0),
    ecas_beh_total_man = conducta_total,
    psychosis_total_man = psicosis,

    # Other (ALS-FTD-Q)
    nps_cognitive_test___99 = if_else(!is.na(alsftdq_total), 1, 0),
    nps_cognitive_test_other = "ALS-FTD-Q",
    nps_cognitive_test_other_score = alsftdq_total,

    # Other fields
    cog_update = today(),
    cognitive_assessment_complete = 2 # Complete
  )

output_samples <- bind_rows(
  output_patient_ids |>
    inner_join(
      biobanco_muestras_ela,
      by = "nhc", relationship = "one-to-many"
    ) |>
    transmute(
      record_id,
      sample_collection_date = coalesce(
        fecha_muestra, fecha_donacion_recepcion, fecha_de_entrada
      ),
      sample_type = tipo_muestra |> case_match(
        "04-Líquido Cefalorraquídeo" ~ 3, # CSF
        "08-Plasma" ~ 2,                  # Serum/Plasma
        "09-Suero" ~ 2,                   # Serum/Plasma
        "13-ADN" ~ 1,                     # DNA
        "Pellet" ~ 1,                     # DNA
        "11-Buffy coat" ~ 1,              # DNA
        "Plasma Litio" ~ 2,               # Serum/Plasma
        "01a-Sangre EDTA" ~ 2,            # Serum/Plasma
        "05-Orina" ~ 4,                   # Urine
        "02-Sangre ocre" ~ 2,             # Serum/Plasma
        .default = 99                     # Other
      ),
      biobank_code = codigo_muestra_noray_banks
    ),
  output_patient_ids |>
    inner_join(
      biobanco_muestras_pals,
      by = "nhc", relationship = "one-to-many"
    ) |>
    transmute(
      record_id,
      sample_collection_date = data_obtencio,
      sample_type = tipus_de_mostra |> case_match(
        c("SÈRUM", "PLASMA") ~ 2, # Serum/Plasma
        "LCR" ~ 3,                # CSF
        .default = 99             # Other
      ),
      biobank_code = codi_mostra
    )
) |>
  mutate(
    sample_storage = 1, # Biobank
    biobank_name = "HUB-ICO-IDIBELL Biobank",
    biobank_link = "https://idibell.cat/en/services/scientific-and-technical-services/biobank/",
    spl_update = today(),
    samples_complete = 2 # Complete
  ) |>
  arrange(sample_collection_date, sample_type, biobank_code)

output_dir <- here("output", str_glue("20260221-{redcap_site_id}"))
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

output_patient_ids |>
  write_csv(
    file.path(output_dir, "redcap-patient-ids.csv"),
    na = ""
  )

list(
  output_status,
  output_consent,
  output_personalinformation,
  output_demographics,
  output_clinicaldata,
  output_diagnosis,
  output_treatment
) |>
  reduce(full_join, by = "record_id") |>
  make_redcap_table(event_name = "basic_arm_1") |>
  filter(!(record_id %in% output_tofersen_study_ids)) |>
  write_redcap_table(file.path(output_dir, "redcap-form-baseline.csv"))

output_alstreatmentdata |>
  make_redcap_table(
    event_name = "visits_arm_1",
    repeat_instrument = "als_treatment_data"
  ) |>
  filter(!(record_id %in% output_tofersen_study_ids)) |>
  write_redcap_table(file.path(output_dir, "redcap-form-alstreatmentdata.csv"))

output_alsfrs |>
  make_redcap_table(
    event_name = "visits_arm_1",
    repeat_instrument = "alsfrs"
  ) |>
  write_redcap_table(file.path(output_dir, "redcap-form-alsfrs.csv"))

output_kingsscale |>
  make_redcap_table(
    event_name = "visits_arm_1",
    repeat_instrument = "kings_scale"
  ) |>
  write_redcap_table(file.path(output_dir, "redcap-form-kingsscale.csv"))

output_eq5d5l |>
  make_redcap_table(
    event_name = "visits_arm_1",
    repeat_instrument = "eq_5d_5l"
  ) |>
  write_redcap_table(file.path(output_dir, "redcap-form-eq5d5l.csv"))

output_fvcinfo |>
  make_redcap_table(
    event_name = "visits_arm_1",
    repeat_instrument = "fvc_info"
  ) |>
  write_redcap_table(file.path(output_dir, "redcap-form-fvcinfo.csv"))

output_advancedrespiratorydata |>
  make_redcap_table(
    event_name = "visits_arm_1",
    repeat_instrument = "advanced_respiratory_data"
  ) |>
  write_redcap_table(file.path(output_dir, "redcap-form-advancedrespiratorydata.csv"))

output_weightandbmi |>
  make_redcap_table(
    event_name = "visits_arm_1",
    repeat_instrument = "weight_and_bmi"
  ) |>
  write_redcap_table(file.path(output_dir, "redcap-form-weightandbmi.csv"))

output_cognitiveassessment |>
  make_redcap_table(
    event_name = "visits_arm_1",
    repeat_instrument = "cognitive_assessment"
  ) |>
  write_redcap_table(file.path(output_dir, "redcap-form-cognitiveassessment.csv"))

output_geneticstudy |>
  make_redcap_table(
    event_name = "samples_analytics_arm_1",
    repeat_instrument = "genetic_study"
  ) |>
  filter(!(record_id %in% output_tofersen_study_ids)) |>
  write_redcap_table(file.path(output_dir, "redcap-form-geneticstudy.csv"))

output_geneticvariations |>
  make_redcap_table(
    event_name = "samples_analytics_arm_1",
    repeat_instrument = "genetic_variations"
  ) |>
  filter(!(record_id %in% output_tofersen_study_ids)) |>
  write_redcap_table(file.path(output_dir, "redcap-form-geneticvariations.csv"))

output_samples |>
  make_redcap_table(
    event_name = "samples_analytics_arm_1",
    repeat_instrument = "samples"
  ) |>
  write_redcap_table(file.path(output_dir, "redcap-form-samples.csv"))
