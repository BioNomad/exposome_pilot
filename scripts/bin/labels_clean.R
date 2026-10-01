label_map <- tibble::tribble(
  ~variable,            ~clean_name,
  # Dust Allergens
  "Bla_g_1_1",          "German Cockroach (dust)",
  "Can_f_1_1",          "Dog (dust)",
  "Der_f_1_1",          "Dust Mite (dust)",
  "Fel_d_1_1",          "Cat (dust)",
  "mus_m_1_1",          "Mouse (dust)",
  
  # IgE — Respiratory
  "sige_cat",           "Cat Dander IgE",
  "sige_cockroach",     "German Cockroach IgE",
  "sige_d_far",         "D. farinae IgE",
  "sige_d_pter",        "D. pteronyssinus IgE",
  "sige_dog",           "Dog IgE",
  "sige_grass",         "Grass Mix IgE",
  "sige_mold",          "Mold Mix IgE",
  "sige_mouse",         "Mouse IgE",
  "sige_ragweed",       "Common Ragweed IgE",
  "sige_tree",          "Tree Mix IgE",
  
  # IgE — Food
  "sige_cow_milk",      "Cow's Milk IgE",
  "sige_egg_white",     "Egg White IgE",
  "sige_peanut",        "Peanut IgE",
  "sige_shrimp",        "Shrimp IgE",
  
  # IgE — Superantigen
  "sige_s_enter_a",     "S. Enterotoxin A IgE",
  "sige_s_enter_b",     "S. Enterotoxin B IgE",
  "sige_s_enter_c",     "S. Enterotoxin C IgE",
  "sige_s_enter_tsst",  "S. Enterotoxin TSST IgE",
  
  # IgE — Total
  "sige_total_ige",     "Total IgE",
  
  # Indoor Air
  "no2f",               "NO2",
  "pm10f",              "PM10",
  "pm25f",              "PM2.5",
  "ufpln",              "Ultrafine Particles",
  "iv_tempavg",         "Indoor Temperature",
  "airnicd",            "Air nicotine",
  
  # Parabens
  "b_pb_ug_g_cr",       "Butyl Paraben",
  "e_pb_ug_g_cr",       "Ethyl Paraben",
  "m_pb_ug_g_cr",       "Methyl Paraben",
  "p_pb_ug_g_cr",       "Propyl Paraben",
  
  # Antimicrobials
  "bp_3_ug_g_cr",       "Benzophenone-3",
  "tcc_ug_g_cr",        "Triclocarban",
  "trcs_ug_g_cr",       "Triclosan",
  
  # Phenols
  "bpa_ug_g_cr",        "Bisphenol A",
  "bpf_ug_g_cr",        "Bisphenol F",
  "bps_ug_g_cr",        "Bisphenol S",
  "x24_dcp_ug_g_cr",    "2,4-Dichlorophenol",
  "x25_dcp_ug_g_cr",    "2,5-Dichlorophenol",
  
  # Serum Essential Metals
  "serum_Cu",           "Copper (serum)",
  "serum_Fe",           "Iron (serum)",
  "serum_Mn",           "Manganese (serum)",
  "serum_Se",           "Selenium (serum)",
  "serum_Zn",           "Zinc (serum)",
  
  # Serum Non-Essential Metals
  "serum_Al",           "Aluminum (serum)",
  "serum_As",           "Arsenic (serum)",
  "serum_Be",           "Beryllium (serum)",
  "serum_Cd",           "Cadmium (serum)",
  "serum_Cr",           "Chromium (serum)",
  "serum_Ni",           "Nickel (serum)",
  "serum_Pb",           "Lead (serum)",
  "serum_Li",           "Lithium (serum)",
  
  # Urine Essential Metals
  "urine_Co",           "Cobalt (urine)",
  "urine_Cu",           "Copper (urine)",
  "urine_Fe",           "Iron (urine)",
  "urine_Mn",           "Manganese (urine)",
  "urine_Mo",           "Molybdenum (urine)",
  "urine_Se",           "Selenium (urine)",
  "urine_Zn",           "Zinc (urine)",
  
  # Urine Non-Essential Metals
  "urine_Al",           "Aluminum (urine)",
  "urine_As",           "Arsenic (urine)",
  "urine_Be",           "Beryllium (urine)",
  "urine_Cd",           "Cadmium (urine)",
  "urine_Cr",           "Chromium (urine)",
  "urine_Li",           "Lithium (urine)",
  "urine_Ni",           "Nickel (urine)",
  "urine_Pb",           "Lead (urine)",
  "urine_Sb",           "Antimony (urine)",
  "urine_V",            "Vanadium (urine)",
)
